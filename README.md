# Parallel-distributed-K-means-algorithm
The idea in this program is to adopt a distributed memory viewpoint of the k-means algorithm https://en.wikipedia.org/wiki/K-means_clustering

## Table of contents
* [General Information](#general-information)
* [Screenshots](#screenshots)
* [Code Examples](#code-examples)
* [Contact](#contact)

## General Information
main.c is the main driver. The datasets used are US pollution data from 2016, which is available on umn cselab machine. Here, I only upload the very small datesets pollution_small.csv for testing. The auxil1.c contains a set of auxilliary functions, which includes the parallel dirstirbuted version of k-means algorithm.

## Screenshots
This is a demo of running
```bash
mpirun -np 8 -hostfile hostfile -map-by node  main.ex
```
on phiXX.cselabs.umn.edu clusters.
![Demo](./img/demo.png)

## Code Examples
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
mpirun -np 16 -hostfile hostfile -map-by node  main.ex
```
where X = 16, <filename> is hosts, and -map-by node will load balance the processes across the available nodes, numbering each process in a round-robin fashion.

To verify the correctness of the program, run:
```bash
python3 solver.py
```
The output should be
```bash
number of correct labels =  1021 / 1021
```
## Contact
Created by [liang664@umn.edu] - feel free to contact me if you have any questions!
