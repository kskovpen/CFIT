#/bin/env bash

PYTHONINC="/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/python/2.7.11-oenich2/include/python2.7/"
BDIR=$(pwd)

swig -c++ -python cfit.i

g++ -fPIC -c -std=c++11 -Wno-deprecated-declarations -m64 cfit_wrap.cxx -I${PYTHONINC} -I$(root-config --cflags) -I${BDIR}

g++ cfit_wrap.o -shared -fPIC -O3 -pthread -m64 -std=c++11 -o _cfit.so -L${BDIR} -lCFIT
