package edu.berkeley.bid.comm;

import java.io.Serializable;
import java.lang.reflect.Array;
import java.util.ArrayList;

/**
 * Created by qualiali on 26/02/2017.
 */
public class WorkerProgress implements Serializable {

    int sentSize = 0;
    int recvSize = 0;
    int lastSentSize = 0;
    int lastRecvSize = 0;

    ArrayList<Long> timestamps = new ArrayList<>();
    ArrayList<Double> sentBandwidthHistory = new ArrayList<>();
    ArrayList<Double> recvBandwidthHistory = new ArrayList<>();

    public WorkerProgress(){
	timestamps.add(System.currentTimeMillis());
    	sentBandwidthHistory.add(0.0);
	recvBandwidthHistory.add(0.0);
    }

    public void addSentBandwidth(double bandwidth){
        sentBandwidthHistory.add(bandwidth);
    }

    public void addRecvBandwidth(double bandwidth){
        recvBandwidthHistory.add(bandwidth);
    }

    public void addSend(int size) {
        sentSize += size;
    }

    public void addRecv(int size){
        recvSize += size;
    }

    public void addTimestamp (long timestamp){
        long lastTimestamp = timestamps.get(timestamps.size()-1);
        timestamps.add(timestamp);
        long interval = timestamp - lastTimestamp;
        int sentSizeIncrement = sentSize - lastSentSize;
        int recvSizeIncrement = recvSize - lastRecvSize;
        lastSentSize = sentSize;
        lastRecvSize = recvSize;
        double sentBandwidth = sentSizeIncrement/1024.0/1024/(interval/1000.0);
        double recvBandwidth = recvSizeIncrement/1024.0/1024/(interval/1000.0);
        sentBandwidthHistory.add(sentBandwidth);
        recvBandwidthHistory.add(recvBandwidth);
    }

    public String toString() {
        String output = "";
        output += "Worker Progress:";

        long lastTimestamp = timestamps.get(timestamps.size()-1);
        long millis = lastTimestamp / 1000000;
        long second = (lastTimestamp / 1000) % 60;
        long minute = (lastTimestamp / (1000 * 60)) % 60;
        long hour = (lastTimestamp / (1000 * 60 * 60)) % 24;
        String time = String.format("%02d:%02d:%02d:%d", hour, minute, second, millis);
        double sentBandwidth = sentBandwidthHistory.get(sentBandwidthHistory.size()-1);
        double recvBandwidth = recvBandwidthHistory.get(recvBandwidthHistory.size()-1);
        output += String.format("Time: %s, Sent: %.2f MB/s, Recv: %.2f Mb/s", time, sentBandwidth, recvBandwidth);
        return output;
    }

    public void printProgress(){
        System.out.println(this.toString());
    }
}
