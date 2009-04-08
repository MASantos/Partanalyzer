			./partanalyze (Partition Analyzer)
                                                                       
./partanalyze aims at being a general program for analyzing (sets of) partitions.
 Here a partition is defined as in set theory of mathematics (see      
 http://en.wikipedia.org/wiki/Partition_of_a_set). It also allows to   
 edit (rudimentarily), as well as generate, partitions.                
                                                                       
Examples:
For lastest options check the help from the program
./partanalyze -h

Check consistency of a given partition test.subfam.lst based on a matrix of interactions given by test-blast_pairwise_id.
How large are the intra-cluster values compared to the inter-cluster ones.
./partanalyze -c test-blast_pairwise_id test.subfam.lst
or
./partanalyze --check-consistency-of-partition test-blast_pairwise_id test.subfam.lst
which also accepts an abreviated form as
./partanalyze --ccop test-blast_pairwise_id test.subfam.lst

Calculate VI distance between two partitions and between each of them and their intersection
 Definition of VI distance: Given two partitions P1 and P2, with cluster size distributions
 {n_k} and {n_k'} respectively, where k and k' are indexes to each of their corresponding clusters, and 
 such that Sum_k n_k = Sum_k n'_k = N, the VI distance is defined as
 
	VI (P1,P2)  = Sum_k n_k/N * log( n_k/N) + Sum_k' n_k'/N * log( n_k'/N) - 2 * Sum_k Sum_k' n_kk'/N log(n_kk'/N)
 where  n_kk' is the number of items common to cluster k of P1 and cluster k' of P2.
 This definition satisfies the triangular inequality, i.e., for any three partitions P1,P2 and P, it is
 	VI (P1, P) + VI (P,P2) >= VI (P1,P2)
./partanalyze --vi-distance test.subfam.lst test.subfam.lst2
or simply
./partanalyze -v test.subfam.lst test.subfam.lst2

Print the intersection of 2 partitions test.subfam.lst and test.subfam.lst2
./partanalyze -i test.subfam.lst test.subfam.lst2
Performs the intersection of P1 and P2 as induced by the intersection operation on
 the underlying set (the one that contains all elements). This gives a new partition
 I such that each cluster of I is obtained as an intersection of one cluster of P1 and
 one of P2 (all againts all).

Print the purity scores for partition1 (target) againts partition2 (reference)
./partanalyze --purity-scores test.subfam.lst test.subfam.lst2
or simplply
./partanalyze -p test.subfam.lst test.subfam.lst2

It outputs the purity strict and purity lax values. 
Purity strict of P1 againts P2 := the number of non-singleton clusters of P1 that
 are exactly identical to one of P2, divided by the number of non-singleton clusters
 of P2 (the reference).
Purity Lax of P1 againts P2 := the number of non-singleton clusters of P1 that
 are subsets of a cluster of P2, divided by the number of non-singleton clusters
 of P1 (the target).

For debugging: print the interaction matrix read by the program
./partanalyze --print-matrix test-blast_pairwise_id


Full help (as of version _VERSION_)
