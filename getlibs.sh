#!/bin/bash

source=$1


BIDMAT_ROOT="${BASH_SOURCE[0]}"
if [ ! `uname` = "Darwin" ]; then
  BIDMAT_ROOT=`readlink -f "${BIDMAT_ROOT}"`
else 
  while [ -L "${BIDMAT_ROOT}" ]; do
    BIDMAT_ROOT=`readlink "${BIDMAT_ROOT}"`
  done
fi
BIDMAT_ROOT=`dirname "$BIDMAT_ROOT"`
BIDMAT_ROOT=`pwd`
BIDMAT_ROOT="$( echo ${BIDMAT_ROOT} | sed s+/cygdrive/c+c:+ )" 

cp ${source}/lib/*.jar ${BIDMAT_ROOT}/lib
cp ${source}/lib/*.so ${BIDMAT_ROOT}/lib
cp ${source}/lib/*.dll ${BIDMAT_ROOT}/lib
cp ${source}/lib/*.dylib ${BIDMAT_ROOT}/lib
cp ${source}/lib/*.jnilib ${BIDMAT_ROOT}/lib

cp ${source}/BIDMat.jar ${BIDMAT_ROOT}

