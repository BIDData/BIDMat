package edu.berkeley.bid;
import java.io.*;
import java.util.zip.*;

public final class UTILS {

    private UTILS() {}

    static {
      LibUtils.loadLibrary("bidmatmkl");
    }
    
    public static native  int  getnumthreads( );
    public static native  void setnumthreads( int n );

    public static native  void memcpybi( int n, byte [] a, long startA, int [] b, long startB );
    public static native  void memcpybf( int n, byte [] a, long startA, float [] b, long startB );
    public static native  void memcpybd( int n, byte [] a, long startA, double [] b, long startB );

    public static native  void memcpyib( int n, int [] a, long startA, byte [] b, long startB );
    public static native  void memcpyfb( int n, float [] a, long startA, byte [] b, long startB );
    public static native  void memcpydb( int n, double [] a, long startA, byte [] b, long startB );
    
    public static native  void memcpyfi( int n, float [] a, long startA, int [] b, long startB );
    public static native  void memcpyif( int n, int [] a, long startA, float [] b, long startB );
    
    public static native  void memcpyli( int n, long [] a, long startA, int [] b, long startB );
    public static native  void memcpyil( int n, int [] a, long startA, long [] b, long startB );
    
    public static native  void lsort(long [] A, int N, int asc);
    
    public static OutputStream _getOutputStream(String fname, int compressionLevel) throws IOException {
    	FileOutputStream fout = new FileOutputStream(fname);
    	switch (compressionLevel) {
    	case 0: return new GZIPOutputStream(fout, 1024*1024){{def.setLevel(0);}};
    	case 1: return new GZIPOutputStream(fout, 1024*1024){{def.setLevel(1);}};
    	case 2: return new GZIPOutputStream(fout, 1024*1024){{def.setLevel(2);}};
    	case 3: return new GZIPOutputStream(fout, 1024*1024){{def.setLevel(3);}};
    	case 4: return new GZIPOutputStream(fout, 1024*1024){{def.setLevel(4);}};
    	case 5: return new GZIPOutputStream(fout, 1024*1024){{def.setLevel(5);}};
    	case 6: return new GZIPOutputStream(fout, 1024*1024){{def.setLevel(6);}};
    	case 7: return new GZIPOutputStream(fout, 1024*1024){{def.setLevel(7);}};
    	case 8: return new GZIPOutputStream(fout, 1024*1024){{def.setLevel(8);}};
    	case 9: return new GZIPOutputStream(fout, 1024*1024){{def.setLevel(9);}};
    	default: throw new RuntimeException("Unsupported compression level "+compressionLevel);
    	}    		
    }
}