

BIDMat is a very fast matric algebra library. Check the latest <b><a href="https://github.com/BIDData/BIDMach/wiki/Benchmarks">benchmarks</a></b>

The github distribution contains source code only. To get the libraries for your platform, run ./getdevlibs.sh from this directory. Then you can run bidmat with ./bidmat. 

You can build the Java/Scala main jar with sbt (included). There are build scripts for both Scala 2.10 and 2.11. Copy the appropriate one into build.sbt, and then do "./sbt package". You can then run bidmat with ./bidmat (you still need to download the libraries some of which are native).

You can also download an executable bundle from <b><a href="http://bid2.berkeley.edu/bid-data-project/download/">here</a></b>. You will need the libs from there in order to build from a git branch. We use a lot of native code which isn't all available from repos, and you will save a lot of time and headaches by grabbing compiled versions.

The main project page is <b><a href="http://bid2.berkeley.edu/bid-data-project/">here</a></b>.

Documentation is <b><a href="https://github.com/BIDData/BIDMat/wiki">here in the wiki</a></b>

BIDMat is a sister project of BIDMach, a machine learning library, which is 
<b><a href="https://github.com/BIDData/BIDMach">also on github</a></b>

