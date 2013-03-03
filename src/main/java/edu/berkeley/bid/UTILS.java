package edu.berkeley.bid;
import java.io.*;
import java.util.zip.*;

public final class UTILS {

    private UTILS() {}

    static {
        System.loadLibrary("bidmatmkl");
    }
    
    public static native  int  getnumthreads( );
    public static native  void setnumthreads( int n );

    public static native  void memcpybi( int n, byte [] a, long startA, int [] b, long startB );
    public static native  void memcpybf( int n, byte [] a, long startA, float [] b, long startB );
    public static native  void memcpybd( int n, byte [] a, long startA, double [] b, long startB );

    public static native  void memcpyib( int n, int [] a, long startA, byte [] b, long startB );
    public static native  void memcpyfb( int n, float [] a, long startA, byte [] b, long startB );
    public static native  void memcpydb( int n, double [] a, long startA, byte [] b, long startB );
    
    public static OutputStream _getOutputStream(String fname, Boolean compressed, int compressionLevel) throws IOException {
    	FileOutputStream fout = new FileOutputStream(fname);
    	if (compressed) {
    		switch (compressionLevel) {
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
    	} else {
    		return new BufferedOutputStream(fout, 1024*1024);
    	}
    }

}