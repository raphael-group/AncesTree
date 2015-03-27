# AncesTree
AnncesTree is an algorithm for clonal tree reconstruction from multi-sample cancer sequencing data.

## Support

For support using AncesTree, please visit the [AncesTree Google Group](https://groups.google.com/forum/#!forum/ancestree).

## Dependencies

AncesTree is written C++. In addition to a recent C++ compiler, it has the following dependencies:

* [CMake](http://www.cmake.org/) (>= 2.8)
* [Boost](www.boost.org) (>= 1.38)
* [LEMON](http://lemon.cs.elte.hu/trac/lemon) graph library (>= 1.3)
* [CPLEX](http://www.ibm.com/developerworks/downloads/ws/ilogcplex/) (>= 12.0)

[Graphviz](http://www.graphviz.org) is required to visualize the resulting DOT files, but this library is not required for compilation.

## Compilation instructions

To compile AncesTree, execute the following commands from the root of the repository:

    mkdir build
    cd build
    cmake ..
    make
    
In case CMake fails to detect either CPLEX or LEMON, run the following command with adjusted paths:

	cmake \
	-DLIBLEMON_ROOT=~/lemon \
	-DCPLEX_INC_DIR=~/ILOG/cplex/include/ \
	-DCPLEX_LIB_DIR=~/ILOG/cplex/lib/x86-64_osx/static_pic \
	-DCONCERT_LIB_DIR=~/ILOG/concert/lib/x86-64_osx/static_pic \
	-DCONCERT_INC_DIR=~/ILOG/concert/include/ ..
	
## Usage instructions

AncesTree can be run as follows:

    ./ancestree ../data/real/RK26.txt -sol RK26.sol -dot RK26.dot

This will save the solution to `RK26.sol` and a GraphViz visualization of the clonal tree and its mixing in `RK26.dot`.

To obtain a PDF of the tree, run the following command:

    dot -Tpdf RK26.dot -o RK26.pdf

### Parameters

### Input format

### Output format

### Benchmarking