# versalignLib

Author: [Tobias Neumann](mailto:tobias.neumann.at@gmail.com)

Programming language: C/C++11

Technologies: [SSE2](https://en.wikipedia.org/wiki/SSE2), [AVX2](https://en.wikipedia.org/wiki/Advanced_Vector_Extensions), [OpenCL](https://www.khronos.org/opencl/), [OpenMP](http://www.openmp.org)

Source: [tar](https://github.com/t-neumann/versalignLib/archive/v0.1.tar.gz), [zip](https://github.com/t-neumann/versalignLib/archive/v0.1.zip) 

What is it?
-----------

**versalignLib** is a showcasing project of parallelization technologies comprising implementations of the sequence alignment algorithms
[Smith-Waterman](https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm)
and [Needleman-Wunsch](https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm).

It builds a set of libraries, each library containing implementations utilizing a different
parallelization technology (OpenMP, SSE2, AVX2, OpenCL). They are provided as [shared objects](https://en.wikipedia.org/wiki/Dynamic_loading)
and can be dynamically loaded depending on available resources and can be benchmarked to one another.

**For any details on usage and benchmarks of the reference implementation, please refer to the [Github wiki](https://github.com/t-neumann/versalignLib/wiki).**

INSTALLATION
============

Building **versalignLib** requires *[cmake](http://www.cmake.org/)* (>=2.8.11), *g++* and an OpenCL implementation.

### Build tools

Typically *cmake* and *g++* should be already available on your machine. If not, please install them.

### OpenMP

To activate OpenMP parallelization, your compiler needs to have OpenMP support. Otherwise all Kernels except the OpenCL Kernel will run single-threaded. OpenMP is natively supported by *gcc* and for newer versions of *Clang*.

### OpenCL

To build the OpenCL library, an OpenCL implementation must be available on your machine.

#### Mac

Mac OS X comes natively with an [OpenCL implementation](https://developer.apple.com/opencl/). No need for further setup.

#### Linux
For Linux, **versalignLib** utilizes the [AMD OpenCL™ APP SDK](http://developer.amd.com/appsdk). All you need to to is to make **versalignLib** aware of the AMD OpenCL™ APP SDK library by setting the following environment variables:

 
```
cd versalignLib
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$PWD/opencl-AMD-sdk/x86_64/lib"
export OPENCL_VENDOR_PATH="$PWD/opencl-AMD-sdk/x86_64/lib/vendors"
```

In case the symlink to the `libOpenCL.so.1` library is not preserved, recreate it like this:

```
(cd opencl-AMD-sdk/x86_64/lib ; rm libOpenCL.so ; ln -s libOpenCL.so.1 libOpenCL.so)
```

### AVX2

*cmake* automatically checks whether AVX2 build support is available and will only produce the shared object if if detects support.

**Caveat 1:** You still need to check for AVX2 instruction support during runtime, for more see the **versalignLib** example implementation!

**Caveat 2:** OpenMP support for the AVX2 Kernel is currently disabled and needs more bugfixing!

### Build

To build **versalignLib** simply run the following commands:

```
cd versalignLib
mkdir build
cd build

# Release
cmake -DCMAKE_BUILD_TYPE=Release ..
# Debug
cmake -DCMAKE_BUILD_TYPE=Debug ..

make
```


System requirements
-------------------

<dl>
<table>
  <tbody>
    <tr>
      <td><b>CPU:</b></td>
      <td>64 bit SSE2 enabled, (optional) AVX2 enabled</td>
    </tr>
    <tr>
      <td><b>RAM:</b></td>
      <td>Tested on regular systems with minimum 4 GB RAM</td>
    </tr>
    <tr>
      <td><b>OS:</b></td>
      <td>Linux (Ubuntu Server 16.04, Debian 7.7), Mac OSX (10.11) </td>
    </tr>
    <tr>
      <td><b>Software:</b></td>
      <td><a href="https://cmake.org/"><i>cmake</i></a> (>=2.8.11)</td>
    </tr>
  </tbody>
</table>
</dl>


