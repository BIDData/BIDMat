package edu.berkeley.bid.comm;

import java.io.*;
import java.net.InetSocketAddress;
import java.net.ServerSocket;
import java.net.Socket;
import java.net.SocketException;
import java.nio.ByteBuffer;
import java.util.ArrayList;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.Hashtable;
import java.util.Enumeration;

public class Machine {
    /* Machine Configuration Variables */
    public final int D;                                               // Depth of the network
    public final int M;                                               // Number of Machines
    public final int imachine;                                        // My identity
    public final Groups groups;                                       // group membership indices
    public final int replicate;                                       // replication factor
    public final InetSocketAddress[] workers;                        // IP/socket
    public final boolean doSim;                                       // Simulation on one machine: send messages directly without sockets
    public int trace = 0;                                             // 0: no trace, 1: high-level, 2: everything
    public int sockBase = 50050;                                      // Socket base address
    public int sockOffset = 2;                                        // Leave space for peer and control ports
    public int sendTimeout = 1000;                                    // in msec
    public int recvTimeout = 1000;
    public int configTimeout = 1000;
    public int reduceTimeout = 1000;
    public boolean useLong = false;
    public int round;
    public int maxk = 0;

    public Vec tov;                                                   // Snapshot of the reduced matrix at each layer going down
    public Vec fromv;                                                 // Snapshot of the result matrix coming up
    public Layer[] layers;                                           // All the layers
    public ByteBuffer[] sendbuf;                                     // buffers, one for each destination in a group
    public ByteBuffer[] recbuf;
    public IVec finalMap;                                             // Map from to --> up at layer D-1
    public Msg[][] messages;                                         // Message queue
    public boolean[][] msgrecvd;                                     // Receiver status
    public boolean[][] amsending;                                    // Sender status
    public ExecutorService executor;
    public ExecutorService sockExecutor;
    public Listener listener;
    public Future<?> listenerFuture;
    public Hashtable<SockWriter, Future<?>> writers;
    public Hashtable<SockReader, Future<?>> readers;
    Network network;

    public WorkerProgress workerProgress = new WorkerProgress();   //Tracks the sock message history

    public Machine(Network p0, Groups groups0, int imachine0, int M0, boolean useLong0, int bufsize, boolean doSim0, int trace0,
                   int replicate0, InetSocketAddress[] workers0) {
        network = p0;
        M = M0;
        doSim = doSim0;
        imachine = imachine0;
        groups = groups0;
        replicate = replicate0;
        useLong = useLong0;
        writers = new Hashtable<SockWriter, Future<?>>();
        readers = new Hashtable<SockReader, Future<?>>();
        if (workers0 == null) {
            workers = new InetSocketAddress[M * replicate];
            for (int i = 0; i < M * replicate; i++)
                workers[i] = new InetSocketAddress("localhost", sockBase + i * sockOffset);
        } else {
            workers = workers0;
        }
        D = groups.depth();
        trace = trace0;
        layers = new Layer[D];
        round = 0;
        int cumk = 1;
        maxk = 1;
        int cumPos = 0;
        for (int level = 0; level < D; level++) {
            int k = groups.nodesInGroup(imachine, level).length;
            maxk = Math.max(maxk, k);
        }
        sendbuf = new ByteBuffer[maxk];
        recbuf = new ByteBuffer[maxk];
        for (int i = 0; i < maxk; i++) {
            sendbuf[i] = ByteBuffer.wrap(new byte[4 * bufsize]);
            recbuf[i] = ByteBuffer.wrap(new byte[4 * bufsize]);
        }
        messages = new Msg[M * replicate][];
        msgrecvd = new boolean[M * replicate][];
        amsending = new boolean[M * replicate][];
        for (int i = 0; i < M * replicate; i++) {
            messages[i] = new Msg[3 * D];
            msgrecvd[i] = new boolean[3 * D];
            amsending[i] = new boolean[3 * D];
        }
        for (int level = 0; level < D; level++) {
            int k = groups.nodesInGroup(imachine, level).length;
            int pos = groups.posInGroup(imachine, level);
            if (useLong) {
                layers[level] = new LongLayer(this, k, cumk, cumPos, pos, level);
            } else {
                layers[level] = new Layer(this, k, cumk, cumPos, pos, level);
            }
            cumPos = cumPos * k + pos;
            cumk *= k;
        }
    }

    public void start(int maxk) {
        executor = Executors.newFixedThreadPool(maxk + 6); // set to 1 for sequential messaging.
        sockExecutor = Executors.newFixedThreadPool(4 + 4 * maxk);
        for (int level = 0; level < D; level++) {
            layers[level].executor = executor;
        }
        listener = new Listener();
        listenerFuture = sockExecutor.submit(listener);
    }

    public void stop() {
        if (sockExecutor != null) {
            listener.stop(true);
            listenerFuture.cancel(true);
            sockExecutor.shutdownNow();
        }
        if (executor != null) executor.shutdownNow();
    }

    public void waitForComms() {
        waitForReads();
        waitForWrites();
    }

    public void waitForReads() {
        Enumeration<SockReader> se = readers.keys();
        while (se.hasMoreElements()) {
            SockReader s = se.nextElement();
            Future<?> fut = null;
            synchronized (readers) {
                fut = readers.get(s);
            }
            while (fut != null && !fut.isDone()) {
                try {
                    Thread.sleep(1);
                } catch (Exception e) {
                }
            }
            synchronized (readers) {
                readers.remove(s);
            }
        }
    }

    public void waitForWrites() {
        Enumeration<SockWriter> se = writers.keys();
        while (se.hasMoreElements()) {
            SockWriter s = se.nextElement();
            Future<?> fut = null;
            synchronized (writers) {
                fut = writers.get(s);
            }
            while (fut != null && !fut.isDone()) {
                try {
                    Thread.sleep(1);
                } catch (Exception e) {
                }
            }
            synchronized (writers) {
                writers.remove(s);
            }
        }
    }

    public void clearFlags() {
        for (int i = 0; i < M * replicate; i++) {
            for (int j = 0; j < 3 * D; j++) {
                if (messages[i][j] != null) {
                    log(String.format("clearing m %d t %d\n", i, j));
                }
                messages[i][j] = null;
                msgrecvd[i][j] = false;
                amsending[i][j] = false;
            }
        }
    }

    public void config(IVec toi, IVec fromi, int round0) {
//		clearFlags();
        round = round0;
        IVec[] outputs = new IVec[2];
        for (int i = 0; i < D; i++) {
            layers[i].config(toi, fromi, outputs, round0);
            toi = outputs[0];
            fromi = outputs[1];
        }
        finalMap = IVec.mapInds(fromi, toi);
    }

    public void config(int[] toi, int[] fromi, int round) {
        config(new IVec(toi), new IVec(fromi), round);
    }

    public void config(LVec toi, LVec fromi, int round0) {
        round = round0;
        LVec[] outputs = new LVec[2];
        for (int i = 0; i < D; i++) {
            ((LongLayer) layers[i]).config(toi, fromi, outputs, round0);
            toi = outputs[0];
            fromi = outputs[1];
        }
        finalMap = LVec.mapInds(fromi, toi);
    }

    public void config(long[] toi, long[] fromi, int round) {
        config(new LVec(toi), new LVec(fromi), round);
    }


    public Vec reduce(Vec tov0, int stride, int round0) {
        round = round0;
        tov = tov0;
        for (int d = 0; d < D; d++) {
            tov = layers[d].reduceDown(tov, stride, round0);
        }
        fromv = tov.mapFrom(finalMap, stride);
        for (int d = D - 1; d >= 0; d--) {
            fromv = layers[d].reduceUp(fromv, stride, round0);
        }
        if (trace > 1)
            log(String.format("Reduce machine %d result nnz %d out of %d\n", imachine, fromv.nnz(), fromv.size()));
        return fromv;
    }

    public float[] reduce(float[] tov, int stride, int round) {
        return reduce(new Vec(tov), stride, round).data;
    }

    public Vec configReduce(IVec toi, IVec fromi, Vec tov0, int stride, int round0) {
//		clearFlags();
        round = round0;
        tov = tov0;
        IVec[] outputs = new IVec[2];
        for (int d = 0; d < D; d++) {
            tov = layers[d].configReduce(toi, fromi, outputs, tov, stride, round0);
            toi = outputs[0];
            fromi = outputs[1];
        }
        finalMap = IVec.mapInds(fromi, toi);
        fromv = tov.mapFrom(finalMap, stride);
        for (int d = D - 1; d >= 0; d--) {
            fromv = layers[d].reduceUp(fromv, stride, round0);
        }
        if (trace > 1)
            log(String.format("ConfigReduce machine %d result nnz %d out of %d\n", imachine, fromv.nnz(), fromv.size()));
        return fromv;
    }

    public float[] configReduce(int[] toi, int[] fromi, float[] tov, int stride, int round) {
        return configReduce(new IVec(toi), new IVec(fromi), new Vec(tov), stride, round).data;
    }

    public Vec configReduce(LVec toi, LVec fromi, Vec tov0, int stride, int round0) {
        round = round0;
        tov = tov0;
        LVec[] outputs = new LVec[2];
        for (int d = 0; d < D; d++) {
            tov = ((LongLayer) layers[d]).configReduce(toi, fromi, outputs, tov, stride, round0);
            toi = outputs[0];
            fromi = outputs[1];
        }
        finalMap = LVec.mapInds(fromi, toi);
        fromv = tov.mapFrom(finalMap, stride);
        for (int d = D - 1; d >= 0; d--) {
            fromv = layers[d].reduceUp(fromv, stride, round0);
        }
        if (trace > 1)
            log(String.format("ConfigReduce machine %d result nnz %d out of %d\n", imachine, fromv.nnz(), fromv.size()));
        return fromv;
    }

    public float[] configReduce(long[] toi, long[] fromi, float[] tov, int stride, int round) {
        return configReduce(new LVec(toi), new LVec(fromi), new Vec(tov), stride, round).data;
    }


    public class SockWriter implements Runnable {
        int dest;
        Msg msg;

        public SockWriter(int dest0, Msg msg0) {
            msg = msg0;
            dest = dest0;
        }

        public void run() {
            Socket socket = null;
//			log(String.format("M %d W %d Running writer %s\n", imachine, dest, this.toString()));
            try {
                long startTime = System.nanoTime();
                socket = new Socket();
                socket.connect(workers[dest], sendTimeout);
                if (socket.isConnected()) {
                    amsending[dest][msg.tag % (3 * D)] = true;
                    DataOutputStream ostr = new DataOutputStream(socket.getOutputStream());
                    ostr.writeInt(msg.size);
                    ostr.writeInt(msg.sender);
                    ostr.writeInt(msg.tag);
                    //TODO: check this part
                    ostr.write(msg.buf, 0, msg.size * 4);
                    //if there is exception when writing, the record will not be logged
                    long endTime = System.nanoTime();
                    workerProgress.addSend(msg.size*4);
                }
            } catch (Exception e) {
                if (trace > 0)
                    log(String.format("Machine %d round %d problem writing socket " + e + "\n", imachine, round));
            } finally {
                try {
                    if (socket != null) socket.close();
                } catch (Exception e) {
                }
                amsending[dest][msg.tag % (3 * D)] = false;
                synchronized (writers) {
//					log(String.format("M %d W %d Removing writer %s\n", imachine, dest, this.toString()));
                    if (writers.containsKey(this)) {
                        writers.remove(this);
                    } else {
//						log(String.format("M %d W %d Couldnt find writer %s\n", imachine, dest, this.toString()));
                    }
                }
            }
        }
    }

    public class SockReader implements Runnable {
        Socket socket = null;
        int round = 0;

        public SockReader(Socket sock, int round0) {
            socket = sock;
            round = round0;
        }

        public void run() {
            if (trace > 2) log(String.format("Machine %d round %d got a packet\n", imachine, round));
            try {
                DataInputStream istr = new DataInputStream(socket.getInputStream());
                long startTime = System.nanoTime();
                int len = istr.readInt();
                int src = istr.readInt();
                int tag = istr.readInt();
                if (trace > 2)
                    log(String.format("Machine %d round %d got packet src %d tag %d\n", imachine, round, src, tag));
                if (tag / (D * 3) < round) {
                    if (trace > 0) log("out of order packet received\n");
                    istr.close();
                } else {
                    int waiting = 0;
                    while (waiting < reduceTimeout && tag / (D * 3) > round) {
                        Thread.sleep(10);
                        waiting += 10;
                    }

                    int tag0 = tag % (3 * D);
                    if (!msgrecvd[src][tag0]) {
                        Msg msg = new Msg(len, src, imachine, tag);
                        //TODO: check this part
                        istr.readFully(msg.buf, 0, len * 4);
                        long endTime = System.nanoTime();
                        workerProgress.addRecv(len*4);
                        synchronized (Machine.this) {
                            if (!msgrecvd[src][tag0]) {
                                messages[src][tag0] = msg;
                                msgrecvd[src][tag0] = true;
                            }
                        }
                        if (trace > 4)
                            log(String.format("m %d r %d from machine %d set tag %d %d\n", imachine, round, src, tag, tag0));
                    } else {
                        if (trace > 0)
                            log(String.format("SocketReader got a duplicate message to %d from %d with tag %d\n", imachine, src, tag));
                    }
                }
            } catch (Exception e) {
                if (trace > 0)
                    log(String.format("Machine %d round %d Problem reading socket " + e.toString() + "\n", imachine, round));
            } finally {
                try {
                    socket.close();
                } catch (IOException e) {
                }
                synchronized (readers) {
                    if (readers.containsKey(this)) readers.remove(this);
                }
            }
        }
    }

    public class Listener implements Runnable {
        boolean stop = false;
        ServerSocket ss = null;

        public Listener() {
            int socknum = workers[imachine].getPort();
            try {
                ss = new ServerSocket(socknum);
            } catch (Exception e) {
                throw new RuntimeException(String.format("Machine couldnt start socket listener on %d ", socknum) + e);
            }
        }

        public void run() {
            while (!stop) {
                try {
                    SockReader scs = new SockReader(ss.accept(), round);
                    Future<?> fut = sockExecutor.submit(scs);
                    synchronized (readers) {
                        readers.put(scs, fut);
                    }
                } catch (SocketException e) {
                    // This is probably due to the server shutting to. Don't do anything.
                } catch (Exception e) {
                    throw new RuntimeException("Socket listener had a problem " + e);
                }
            }
        }

        public boolean stillSending() {
            boolean sending = false;
            for (int i = 0; i < amsending.length; i++) {
                boolean[] sendrow = amsending[i];
                for (int j = 0; j < sendrow.length; j++) {
                    if (amsending[i][j]) sending = true;
                }
            }
            return sending;
        }

        public void stop(boolean force) {
            if (!force) {
                while (stillSending()) {
                    try {
                        Thread.sleep(1);
                    } catch (InterruptedException e) {
                    }
                }
            }
            try {
                stop = true;
                ss.close();
            } catch (Exception e) {
                throw new RuntimeException("Trouble closing listener");
            }
        }
    }

    public void dumptags(String s) {
        synchronized (network) {
            System.out.print(s);
            for (int i = 0; i < M; i++) {
                for (int j = 0; j < M; j++) {
                    for (int k = 0; k < 3 * D; k++) {
                        if (network.machines[i].messages[j][k] != null) {
                            System.out.print("1");
                        } else {
                            System.out.print("0");
                        }
                    }
                }
            }
            System.out.print("\n");
        }
    }

    public boolean sendrecv(int ith, ByteBuffer sbuf, int sendn, int outi, ByteBuffer rbuf, int recn, int ini, int tag0) {
        if (trace > 4)
            log(String.format("m %d th %d r %d sendrecv enter out %d in %d tag %d\n", imachine, ith, round, outi, ini, tag0));
        if (trace > 2)
            log(String.format("Round %d sendrecv machine %d to %d from %d tag %d\n", round, imachine, outi, ini, tag0));
        sbuf.rewind();
        Msg msg = new Msg(sbuf.array(), sendn, imachine, outi, tag0);
        int tag = tag0 % (3 * D);
        if (imachine == outi) {
            rbuf.clear();
            rbuf.put(msg.buf, 0, 4 * sendn);
            rbuf.rewind();
            if (trace > 2)
                log(String.format("Round %d sendrecv machine %d to %d from %d tag %d done\n", round, imachine, outi, ini, tag0));
            if (trace > 4)
                log(String.format("m %d th %d r %d sendrecv exit out %d in %d tag %d\n", imachine, ith, round, outi, ini, tag0));
            return true;
        } else {
            if (trace > 4)
                dumptags(String.format("m %d th %d r %d enter  %d %d %d ", imachine, ith, round, outi, ini, tag0));
            if (doSim) {
                for (int i = 0; i < replicate; i++) {
                    if (trace > 0) {
                        if (network.machines[outi + i * M].messages[imachine][tag] != null)
                            log(String.format("Round %d sendrecv machine %d to %d from %d tag %d msg exists\n", round, imachine, outi, ini, tag0));
                    }
                    if (trace > 4)
                        log(String.format("m %d th %d r %d to machine %d set tag %d %d\n", imachine, ith, round, outi, tag0, tag));
                    synchronized (this) {
                        network.machines[outi + i * M].messages[imachine][tag] = msg;
                    }
                }
            } else {
                for (int i = 0; i < replicate; i++) {
                    SockWriter w = new SockWriter(outi + i * M, msg);
                    Future<?> fut = sockExecutor.submit(w);
                    synchronized (writers) {
                        writers.put(w, fut);
//						log(String.format("M %d W %d Starting writer %s\n", imachine, outi, fut.toString()));
                    }
                }
            }
            if (trace > 4)
                dumptags(String.format("m %d th %d r %d change %d %d %d ", imachine, ith, round, outi, ini, tag0));
            boolean gotit = false;
            int waiting = 0;
            while (!gotit && waiting < recvTimeout) {
                for (int i = 0; i < replicate; i++) {
                    if (messages[ini + i * M][tag] != null) {
                        Msg rmsg = messages[ini + i * M][tag];
                        rbuf.clear();
                        rbuf.put(rmsg.buf, 0, 4 * rmsg.size);
                        rbuf.rewind();
                        gotit = true;
                        break;
                    }
                    try {
                        Thread.sleep(10);
                        waiting += 10;
                    } catch (InterruptedException e) {
                    }
                }
            }
            Boolean success = true;
            if (trace > 4)
                dumptags(String.format("m %d th %d r %d later  %d %d %d ", imachine, ith, round, outi, ini, tag0));
            if (waiting >= recvTimeout) {
                log(String.format("m %d th %d r %d ", imachine, ith, round));
                log(String.format("Round %d sendrecv machine %d to %d from %d tag %d timed out\n", round, imachine, outi, ini, tag0));
                success = false;
            }
            for (int i = 0; i < replicate; i++) {
                if (trace > 4)
                    log(String.format("m %d th %d r %d from machine %d clear tag %d %d\n", imachine, ith, round, ini, tag0, tag));
                synchronized (this) {
                    messages[ini + i * M][tag] = null;
                    msgrecvd[ini + i * M][tag] = false;
                }
            }
            if (trace > 2)
                log(String.format("Round %d sendrecv machine %d to %d from %d tag %d done after %d\n", round, imachine, outi, ini, tag0, waiting));
            if (trace > 4)
                dumptags(String.format("m %d th %d r %d exit   %d %d %d ", imachine, ith, round, outi, ini, tag0));
            if (trace > 4)
                log(String.format("m %d th %d r %d sendrecv exit out %d in %d tag %d\n", imachine, ith, round, outi, ini, tag0));
            return success;
        }
    }

    public void log(String msg) {
        if (network != null) {
            synchronized (network) {
                System.out.println(msg);
            }
        } else {
                System.out.println(msg);
        }
    }
}