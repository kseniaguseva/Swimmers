CXXFLAGS=-Wall -O99 -g -fPIC  -lpython3.11 -lboost_python3  -lboost_python3 -lboost_numpy3 -I. -I/usr/include/python3.11 -I/usr/lib/python3.11/site-packages/numpy/core/include/  -shared -Wl,-soname,lifluid.so

ALL: libfluid.so

libfluid.so: fluid.cc demangle.cc
	g++ fluid.cc demangle.cc -o libfluid.so ${CXXFLAGS} 
