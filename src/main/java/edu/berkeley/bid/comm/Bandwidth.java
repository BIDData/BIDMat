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
    public long startTime;           //start time in milli seconds
    public long endTime;             //end time in milli seconds

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

        //logging
        System.out.printf("iMachine: %d, round: %d, src: %d, dst: %d, size: %d, startTime: %d, endTime: %d",
                iMachine, round, src, dst, size, startTime, endTime);
    }

}
