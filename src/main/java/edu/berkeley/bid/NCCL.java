package edu.berkeley.bid;
import jcuda.Pointer;
import jcuda.runtime.cudaStream_t;

public final class NCCL {

    private long handle = 0;

    public NCCL() {}

    protected void finalize() {
        if (handle != 0) {
            ncclCommDestroy(this);
            handle = 0;
        }
    }

    static {
//      LibUtils.loadLibrary("bidmatcpu", true);
    }

    public static native int hasNCCL();

    public static native int ncclCommInitAll(NCCL [] comms);

    public static native void ncclCommDestroy(NCCL comm);

    /**
     * Should pass float pointers for send/recv buffers. 
     * Enumeration types are defined here: https://github.com/NVIDIA/nccl/blob/master/src/nccl.h
     */

    public static native int ncclAllReduce(Pointer send, Pointer recv, int count, int datatype, int op, NCCL comm, cudaStream_t stream);
    
}
