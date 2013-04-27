Source code for parallel implementation of Regularized Markov Clustering (RMCL) based on message passing interface. 

Input should be an edge file with src_id"\t"dst_id per line. Vertex ids are assumed to be numbered from 0.
Output file (with name as provided through command line argument) will have cluster_id per line with the following interpretation:
value x on line y means vertex with id y belongs to cluster x.

To compile and run:
make			(Makefile uses mpicxx but mpiCC may also work based on mpi library being used)
./rmcl_mpi --help for command line arguments
mpiexec -n # ./rmcl_mpi	COMMAND LINE ARGUMENTS
mpirun -np # ./rmcl_mpi COMMAND LINE ARGUMENTS

Comments,bugs and questions are welcomed at varia.siddharth@gmail.com.



Note: compiled and tested using openmpi on ubuntu 11.04
