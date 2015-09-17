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

if [ `uname` = "Darwin" ]; then
    binnames=".dylib,.jnilib"
    cdir="osx"
elif [ "$OS" = "Windows_NT" ]; then
    binnames=".dll"
    cdir="win"
else 
    binnames=".so"
    cdir="linux"
fi

source="http://bid2.berkeley.edu/bid-data-project"

cd ${BIDMAT_ROOT}/lib
wget -r -A.txt,.html,.jar ${source}/lib/
wget -r -A${binnames} ${source}/lib/

mv ${BIDMAT_ROOT}/lib/BIDMat.jar ${BIDMAT_ROOT}

