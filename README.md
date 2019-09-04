# xkt
Application to test kernel iteration over different mesh accessors. Application simulates creating an arbitrarily refined mesh in 1,2, or 3D, then runs a kernel (a daxpy) over three different mesh accessors: all cells, most refined cells, and contiguous chunks of refined cells. Application main is written in Fortran, and kernels in C++ to test language interoperability.
