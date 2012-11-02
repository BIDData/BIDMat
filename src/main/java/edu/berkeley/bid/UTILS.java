package edu.berkeley.bid;

public final class UTILS {

    private UTILS() {}

    static {
        System.loadLibrary("bidmatmkl");
    }

    public static native  void memcpybi( int n, byte [] a, int startA, int [] b, int startB );
    public static native  void memcpybf( int n, byte [] a, int startA, float [] b, int startB );
    public static native  void memcpybd( int n, byte [] a, int startA, double [] b, int startB );

    public static native  void memcpyib( int n, int [] a, int startA, byte [] b, int startB );
    public static native  void memcpyfb( int n, float [] a, int startA, byte [] b, int startB );
    public static native  void memcpydb( int n, double [] a, int startA, byte [] b, int startB );

}