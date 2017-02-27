package edu.berkeley.bid.comm;

/**
 * Created by qualiali on 26/02/2017.
 */
public class Bandwidth {
    public int iMachine;           //machine's identity
    public int round;              //
    public int src;                 //src machine no
    public int dst;                 //dst machine no
    public int size;                //size of socket message in bytes
    public long startTime;           //start time in nano seconds
    public long endTime;             //end time in nano seconds
    public long duration;           //time spent on sending the socket, in ns
    public double bandwidth;        //bandwidth, in KB/s

    public Bandwidth(int iMachine,
                     int round,
                     int src,
                     int dst,
                     int size,
                     long startTime,
                     long endTime){
        this.iMachine = iMachine;
        this.round = round;
        this.src = src;
        this.dst = dst;
        this.size = size;
        this.startTime = startTime;
        this.endTime = endTime;
        this.duration = endTime-startTime;
        this.bandwidth = size/1024.0/duration*1000000000;

        //logging
//        System.out.printf("iMachine: %d, round: %d, src: %d, dst: %d, \n" +
//                        "size: %d B, startTime: %d ms, endTime: %d ms\n" +
//                        "duration: %d ms, bandwidth: %f KB/s",
//                iMachine, round, src, dst, size, startTime, endTime, duration, bandwidth);
        System.out.printf("src: %d, dst: %d, size: %d B, duration: %d ms, bandwidth: %f KB/s",
                src, dst, size, duration, bandwidth);
    }

}
