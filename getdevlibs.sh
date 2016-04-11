#!/bin/bash

if [ "${ARCH}" = "" ];then
    ARCH=`arch`
fi

BIDMAT_ROOT="${BASH_SOURCE[0]}"
if [ ! `uname` = "Darwin" ]; then
  BIDMAT_ROOT=`readlink -f "${BIDMAT_ROOT}"`
else 
  while [ -L "${BIDMAT_ROOT}" ]; do
    BIDMAT_ROOT=`readlink "${BIDMAT_ROOT}"`
  done
fi
BIDMAT_ROOT=`dirname "$BIDMAT_ROOT"`
pushd "${BIDMAT_ROOT}"  > /dev/null
BIDMAT_ROOT=`pwd`
BIDMAT_ROOT="$( echo ${BIDMAT_ROOT} | sed s+/cygdrive/c+c:+ )" 

source="http://www.cs.berkeley.edu/~jfc/biddata"
cd ${BIDMAT_ROOT}/lib

if [ `uname` = "Darwin" ]; then
    subdir="osx"
    curl -o liblist.txt ${source}/lib/liblist_osx.txt 
elif [ "$OS" = "Windows_NT" ]; then
    subdir="win"
    curl -o liblist.txt ${source}/lib/liblist_win.txt
else
    if [[ "${ARCH}" == arm* || "${ARCH}" == aarch* ]]; then
        subdir="linux_arm"
        curl -o liblist.txt ${source}/lib/liblist_linux_arm.txt
    else
        subdir="linux"
        curl -o liblist.txt ${source}/lib/liblist_linux.txt
    fi
fi
curl -o exelist.txt ${source}/lib/exelist.txt

while read fname; do
    echo -e "\nDownloading ${fname}"
    curl --retry 2  -z ${fname} -o ${fname} ${source}/lib/${fname}
    chmod 755 ${fname}
done < liblist.txt

mv ${BIDMAT_ROOT}/lib/BIDMat.jar ${BIDMAT_ROOT}
rm ${BIDMAT_ROOT}/lib/BIDMach.jar

cp ${BIDMAT_ROOT}/lib/*bidmat*.{so,dll,dylib} ${BIDMAT_ROOT}/src/main/resources/lib
cp ${BIDMAT_ROOT}/lib/*iomp5*.{so,dll,dylib} ${BIDMAT_ROOT}/src/main/resources/lib

cd ${BIDMAT_ROOT}/src/main/resources
libs=`echo lib/*`

cd ${BIDMAT_ROOT}
jar uvf BIDMat.jar -C src/main/resources $libs

