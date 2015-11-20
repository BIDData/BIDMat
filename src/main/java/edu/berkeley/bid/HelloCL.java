package edu.berkeley.bid;

public final class HelloCL {

  static {
    System.out.println("Loading hello_cl library");
    System.loadLibrary("hello_cl");
    System.out.println("Loaded library?");
  }

  public native void foo();

}
