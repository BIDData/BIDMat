@ECHO OFF
:: Set JAVA_HOME here if not set in environment
:: SET JAVA_HOME= 
:: Set as much memory as possible
(SET JAVA_OPTS=-Xmx12G -Xms128M)
:: Fix these if needed
SET JCUDA_VERSION=0.5.0
SET JCUDA_LIBDIR=%CD%\lib
SET LIBDIR=%CD%\lib
SET PATH=%LIBDIR%;%PATH%


SET BIDMAT_LIBS=%CD%\BIDMat.jar;%LIBDIR%\ptplot.jar;%LIBDIR%\ptplotapplication.jar;%LIBDIR%\jhdf5.jar;%LIBDIR%\commons-math3-3.1.1.jar;%LIBDIR%\lz4-1.1.2.jar

SET JCUDA_LIBS=%JCUDA_LIBDIR%\jcuda-%JCUDA_VERSION%.jar;%JCUDA_LIBDIR%\jcublas-%JCUDA_VERSION%.jar;%JCUDA_LIBDIR%\jcufft-%JCUDA_VERSION%.jar;%JCUDA_LIBDIR%\jcurand-%JCUDA_VERSION%.jar;%JCUDA_LIBDIR%\jcusparse-%JCUDA_VERSION%.jar

SET ALL_LIBS=%BIDMAT_LIBS%;%JCUDA_LIBS%;%JAVA_HOME%\lib\tools.jar
:: echo %ALL_LIBS%

scala -nobootcp -cp "%ALL_LIBS%" -Yrepl-sync -i %LIBDIR%\bidmat_init.scala