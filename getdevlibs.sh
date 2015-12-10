#!/bin/bash

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
    subdir="linux"
    curl -o liblist.txt ${source}/lib/liblist_linux.txt
fi
curl -o exelist.txt ${source}/lib/exelist.txt

while read fname; do
    echo -e "\nDownloading ${fname}"
    curl --retry 2 -O ${source}/lib/${fname}
    chmod 755 ${fname}
done < liblist.txt

mv ${BIDMAT_ROOT}/lib/BIDMat.jar ${BIDMAT_ROOT}
rm ${BIDMAT_ROOT}/lib/BIDMach.jar 

