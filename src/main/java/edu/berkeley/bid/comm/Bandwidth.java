package edu.berkeley.bid.comm;

import java.lang.reflect.Array;
import java.util.ArrayList;

/**
 * Created by qualiali on 26/02/2017.
 */
public class Bandwidth {

    public int totalSize = 0;

    public ArrayList<Record> records = new ArrayList<>();

    class Record{
        public int iMachine;           //machine's identity
        public int round;              //
        public int src;                 //src machine no
        public int dst;                 //dst machine no
        public int size;                //size of socket message in bytes
        public long startTime;           //start time in nano seconds
        public long endTime;             //end time in nano seconds
        public long duration;           //time spent on sending the socket, in ns
        public double bandwidth;        //bandwidth, in KB/s

        public Record(int iMachine,
                     int round,
                     int src,
                     int dst,
                     int size,
                     long startTime,
                     long endTime) {
            this.iMachine = iMachine;
            this.round = round;
            this.src = src;
            this.dst = dst;
            this.size = size;
            this.startTime = startTime;
            this.endTime = endTime;
            this.duration = endTime - startTime;
            this.bandwidth = size / 1024.0 / duration * 1000000000;

            // System.out.printf("src: %d, dst: %d, size: %d B, duration: %d ms, bandwidth: %f KB/s\n",
            //         src, dst, size, duration, bandwidth);
        }
    }

    public void addRecord(int iMachine,
                     int round,
                     int src,
                     int dst,
                     int size,
                     long startTime,
                     long endTime){
        records.add(
            new Record(iMachine, round, src, dst, size, startTime, endTime));
        totalSize += size;
    }

    public void printRecords(){
        for(Record r: records){
            long millis = r.startTime / 1000000;
            long second = (millis / 1000) % 60;
            long minute = (millis / (1000 * 60)) % 60;
            long hour = (millis / (1000 * 60 * 60)) % 24;

            String time = String.format("%02d:%02d:%02d:%d", hour, minute, second, millis);
            System.out.println("time: "+time+" total bytes: "+r.size);
        }
    }
}
