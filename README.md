# versalignLib

Author: [Tobias Neumann](mailto:tobias.neumann.at@gmail.com)

Programming language: C/C++11

Technologies: [SSE2](https://en.wikipedia.org/wiki/SSE2), [AVX2](https://en.wikipedia.org/wiki/Advanced_Vector_Extensions), [OpenCL](https://www.khronos.org/opencl/), [OpenMP](www.openmp.org)

What is it?
-----------

**versalignLib** is a set of libraries of implementations of the sequence alignment algorithms
[Smith-Waterman](https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm)
and [Needleman-Wunsch](https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm).

Each library contains implementations utilizing a different
parallelization technology (OpenMP, SSE2, AVX2, OpenCL). They are provided as [shared objects](https://en.wikipedia.org/wiki/Dynamic_loading)
and can be dynamically loaded depending on available resources and can be benchmarked to one another.

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
      <td><a href="https://cmake.org/"><i>cmake</i></a></td>
    </tr>
  </tbody>
</table>
</dl>


