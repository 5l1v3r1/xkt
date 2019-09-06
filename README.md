# xkt
Application to test kernel iteration over different mesh accessors.
- simulates creating an arbitrarily refined mesh in 1,2, or 3D
- then runs a kernel (a daxpy) over three different mesh accessors
    - all cells
    - most refined cells
    - contiguous chunks of refined cells

The application main is written in Fortran, and kernels in C++ to test language interoperability.

# Getting the Code
To check out the source code for xkt,

```bash
    git clone --recursive https://github.com/lanl/xkt.git
```

## Branches
Consider `master` to be the stable branch. The branch `develop` will contain integrated
features as they move towards release. Then, feature branches as necessary.

# Requirements
- Compiler suite with both Fortran (F08) and C++ (C++11) support.
- A build of [Kokkos](https://github.com/kokkos/kokkos) built with the compiler suite above.
- Up-to-date build of [FLCL](https://github.com/kokkos/kokkos-fortran-interop) (included as a submodule).

# Build Instructions
First, build the fortran/Kokkos interop parts:
- export KOKKOS_ROOT=/path/to/kokkos/installation
- cd $(XTKROOT)/lib/flcl/build
- make libflcl.a
- (optionally, run the tests)
    - make test_flcl.x
    - ./test_flcl.x

Next, build and run xkt itself:
- cd $(XKTROOT)/build
- ./xkt.x


# Usage
xkt is meant to simulate running a kernel of an AMR-like mesh with three differing
access patterns. To wit, it has a few compile-time parameters to control the rough
size of the mesh `m`, the dimension of the mesh `dimension`, and the amount of
refinement, both `percentage` and `refine_steps`.

Sensible values for `m` might be from low tens of thousands to several million.

`dimension` is meant to correspond to a physical dimension, and so should be 1, 2, or 3.

`percentage` is roughly how much of the mesh is refined, and `refine_steps` is how
many times the mesh is refined. All of the refinement is at random locations away from
the boundaries of the mesh.

Since the kernels are all Kokkos `parallel_for`, we will typically use Kokkos'
profiling hooks to monitor performance. See: [Kokkos profiling](https://github.com/kokkos/kokkos-tools)

# Feedback
Please raise an issue using the [GitHub issues](https://github.com/lanl/xkt/issues)
feature.

# Release

C19050 - xkt has been acknowledged by NNSA for open source release.

# Copyright

Copyright &copy; 2019. Triad National Security, LLC. All rights reserved.
 
This program was produced under U.S. Government contract 89233218CNA000001 for
Los Alamos National Laboratory (LANL), which is operated by Triad National
Security, LLC for the U.S. Department of Energy/National Nuclear Security
Administration. All rights in the program are reserved by Triad National
Security, LLC, and the U.S. Department of Energy/National Nuclear Security
Administration. The Government is granted for itself and others acting on
its behalf a nonexclusive, paid-up, irrevocable worldwide license in this
material to reproduce, prepare derivative works, distribute copies to the
public, perform publicly and display publicly, and to permit others to do so.
 
This program is open source under the BSD-3 License.
 
Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:
 
1. Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the
   names of its contributors may be used to endorse or promote products
   derived from this software without specific prior written permission.
 
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.