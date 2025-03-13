CXXFLAGS=-Wall -O3 -g -std=c++17 -fPIC -lpython3.13  -lboost_python3  -lboost_numpy3 -I. -I/usr/include/python3.13 -I/usr/lib/python3.13/site-packages/numpy/_core/include/ -shared -Wl,-soname,libfluid.so

ALL: libfluid.so

libfluid.so: demangle.cc numpy_bind.hh
	g++ fluid.cc demangle.cc -o libfluid.so ${CXXFLAGS} 
