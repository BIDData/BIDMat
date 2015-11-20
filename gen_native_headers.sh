mkdir classes
mkdir headers

CLASSPATH=./classes:$CLASSPATH

javac -d ./classes src/main/java/edu/berkeley/bid/LibUtils.java
javac -d ./classes src/main/java/edu/berkeley/bid/HelloCL.java
javah -d ./headers edu.berkeley.bid.HelloCL
