
BIDMat is a matrix library intended to support large-scale exploratory
data analysis and to accelerate production deployment on single
machines or clusters.  While there are many excellent tools exist to
support data analysis at small scale, there is a dearth of tools that
support large-scale analysis or scale-up. Specifically, some goals of
BIDMat are:

1. To provide an interactive data analysis environment, similar to R or
Matlab. However, since we use the Scala language we have the advantage
of a high-end programming language including good general-purpose data
sructures. And also of Scala's compiler-based REPL (Read-Eval-Print Loop). 

2. To leverage native machine performance through native libraries
(Intel MKL, HDF5 and string/XML processing). Java/Scala are excellent
high-level languages, but are one or two orders of magnitude away from
native performance in some key areas: especially matrix algebra and
string processing, and below the bar to a lesser degree in File-IO.

3. To leverage GPU hardware and GPU-based data as a first-class
object.  GPUs now offer large improvements (again one or more orders
of magnitude) over CPU performance in many areas that are relevant to
data mining: matrix algebra, transcendental functions, random number
generation, network and graph algorithms and even natural language
parsing. Our own work suggests that the list is going to continue to
grow, and that GPU acceleration will fairly soon be a requirement for
competitive performance in most algorithms.
 
