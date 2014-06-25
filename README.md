BIDMat is a matrix library intended to support large-scale exploratory
data analysis. Its sister library <b><a href="https://github.com/BIDData/BIDMach">BIDMach</a></b>
implements the machine learning layer.
BIDMat is designed to be both fast and usable. Its
speed comes from a combination of hardware-accelerated (CPU and GPU)
kernels and roofline design: kernels have been systematically
optimized against theoretical limits.  See the <b><a href="https://github.com/BIDData/BIDMach/wiki/Benchmarks">Benchmarks Page</a></b>. BIDMat is written in the Scala
language which is ideal for an interactive, math-like DSL.
The goals of BIDMat are:

1. To provide an interactive data analysis environment, similar to R
or scipy or Matlab. Thanks to the Scala language, we inherit a
REPL, a very flexible operator syntax (including unicode math
operators), and a very clean concurrent programming model using
actors.

2. To leverage native machine performance through native libraries
(Intel MKL, HDF5, CUDA and string/XML processing). Java/Scala are
excellent high-level languages, but are one or two orders of magnitude
away from native performance in some key areas: especially matrix
algebra and string processing, and below the bar to a lesser degree in
File-IO.

3. To leverage GPU hardware and GPU-based data as a first-class
object.  GPUs now offer large improvements (again one or more orders
of magnitude) over CPU performance in many areas that are relevant to
data mining: matrix algebra, transcendental functions, random number
generation. These advantages in low-level operations carry over to
network and graph algorithms and even natural language parsing. 

<b><a href="https://github.com/BIDData/BIDMat/wiki">Documentation is here</a></b>
on the github wiki page.
 
