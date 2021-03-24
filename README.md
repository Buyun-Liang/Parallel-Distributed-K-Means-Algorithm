## Parallel-distributed-K-means-algorithm
The idea in this program is to adopt a distributed memory viewpoint of the k-means algorithm https://en.wikipedia.org/wiki/K-means_clustering

## Files Description
main.c is the main driver. The datasets used are US pollution data from 2016, which is available on umn cselab machine. auxil1.c contains a set of auxilliary functions, which includes the parallel dirstirbuted version of k-means algorithm.

## Setup
To build the executable program main.ex, run:
```bash
make
```

The program could be run on the phiXX.cselabs.umn.edu cluster by running:
```bash
mpirun [ -np X ] [ --hostfile <filename> ] <program>
```
where X is the number of processors, and -hostfile is used for specifying host nodes. For example:
```bash
mpirun -np 16 -hostfile hosts -map-by node  main.ex
```
where X = 16, <filename> is hosts, and -map-by node will load balance the processes across the available nodes, numbering each process in a round-robin fashion.