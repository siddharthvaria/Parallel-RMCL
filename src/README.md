Source code for parallel implementation of Regularized Markov Clustering (RMCL) based on message passing interface. 

RMCL is a graph clustering algorithm based on flow simulations. It works for undirected graphs only. For more info about the algorithm contact me.

Input should be an edge file with src_id"\t"dst_id per line. Each edge should appear only once in the input file, so for an edge from 3 to 5,input file should either 3\t5 or 5\t3 and not both. 

Vertex ids are assumed to be numbered from 0.

Output file (with name as provided through command line argument) will have cluster_id per line with the following interpretation:

value x on line y means vertex with id y belongs to cluster x.

To compile and run:

make			(Makefile uses mpicxx but mpiCC may also work based on mpi library being used)

./rmcl_mpi --help for command line arguments

mpiexec -n # ./rmcl_mpi	COMMAND LINE ARGUMENTS

mpirun -np # ./rmcl_mpi COMMAND LINE ARGUMENTS


Comments,bugs and questions are welcomed at varia.siddharth@gmail.com.



Note: compiled and tested using openmpi on ubuntu 11.04
