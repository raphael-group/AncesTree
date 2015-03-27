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

[Graphviz](http://www.graphviz.org) is required to visualize the resulting DOT files, but is not required for compilation.

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

The `ancestree` exectuable takes the following arguments as input:

	./ancestree [--alpha|-a num] [--beta|-b num] [--dot|-d str]
	   [--gamma|-g num] [--help|-h|-help] [--sol|-s str] [--time|-t int]
	   [--version|-v] read_count_file
	   
where

ARGUMENT       | DEFAULT | DESCRIPTION                                                
---------------|---------|-------------------------------------------------------------
--alpha/-a     | 0.3     | Controls the clustering of mutations in the graph clustering phase: only arcs (v_j, v_k) with 0.5 - alpha <= min_p P(X_pj < X_pk) <= 0.5 + alpha  are considered
--beta/-b      | 0.8     | Controls the confidence in ancestral relationships in the graph: there is an arc (v_j, v_k) if min_p P(X_pj < X_pk) >= beta
--gamma/-g     | 0.01    | Controls the allowed pertubation of observed variant frequencies by defining (1 - gamma) confidence intervals 
--dot/-d       |         | DOT output filename for the clonal tree visualization
--sol/-s       | STDOUT  | Solution output filename
--time/-t      | -1      | ILP time limit in seconds, use -1 for no time limit
--help/-h      |         | Shows usage instructions
--version/-v   |         | Shows version number
read_count_file|         | Input file containing read counts

### Example

To run AncesTree on patient RK26, do:

    ./ancestree ../data/real/CLL077_whole.txt --sol CLL077_whole.sol --dot CLL077_whole.dot

This will save the solution to [CLL077_whole.sol](doc/CLL077_whole.sol) and a Graphviz visualization of the clonal tree and its mixing to [CLL077_whole.dot](doc/CLL077_whole.dot). See below for details on the [input](#input-format) and [output](#output-format) format. 

To obtain a PNG of the tree, run the following command:

    dot -Tpng CLL077_whole.dot -o CLL077_whole.png

![CLL077_whole.png](doc/CLL077_whole.png)

### Input format

The input is a tab-separated ASCII text file. The first line contains the sample headers. The first column contains gene ids. Then every consecutive pair of columns contains read counts for reference alleles and alternate alleles, respectively.

	gene_id	a	a	b	b	c	c	d	d	e	e
	C3orf43	16	13	28	17	35	24	21	22	30	33
	CNOT7	29	17	27	22	21	24	25	22	15	24
	IRF4	36	4	30	10	33	8	25	11	22	13

### Output format



### Benchmarking