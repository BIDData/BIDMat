

BIDMat is a very fast matric algebra library. Check the latest <b><a href="https://github.com/BIDData/BIDMach/wiki/Benchmarks">benchmarks</a></b>

The github distribution contains source code only. To build the system, you need a Java JDK 8, an installation of CUDA 8.0 (if you want to use NVIDIA GPUs), and a copy of <b><a href="https://maven.apache.org/download.cgi">apache maven 3.x</a></b>. On windows, you also need a unix command package like <b><a href="https://www.cygwin.com/">cygwin</a></b>. With those prerequisites, you can do:
<pre>
mvn clean install
</pre>
to build and install and then
<pre>
./bidmat
</pre>
To start bidmat. 

The main project page is <b><a href="http://bid2.berkeley.edu/bid-data-project/">here</a></b>.

Documentation is <b><a href="https://github.com/BIDData/BIDMat/wiki">here in the wiki</a></b>

BIDMat is a sister project of BIDMach, a machine learning library, which is 
<b><a href="https://github.com/BIDData/BIDMach">also on github</a></b>


Take a look at BIDMach_RL, a new project on reinforcement learning which has state-of-the-art implementations of several RL algorithsms: <b><a href="https://github.com/BIDData/BIDMach">also on github</a></b>

