# Maximum Expected Accuracy Partial Order Alignment

Java implementation of a Maximum Expected Accuracy Partial Order Alignment method.

This is the core alignment module that I worked on during my PhD - it takes a set of unaligned sequences and aligns them 
together using a pair-HMM model in a way similar to methods implemented in [PRANK](http://wasabiapp.org/software/prank/ "PRANK")
 and [PAGAN](http://wasabiapp.org/software/pagan/phylogenetic_multiple_alignment/ "PAGAN").

The difference in my implementation is that we don't use a Viterbi algorithm for computing the path through an alignment
and the order in which sequences should be aligned. Instead, we use a Maximum Expected Accuracy algorithm, similar to the 
algorithm employed in consistency-based aligners such as [PROBCONS](http://probcons.stanford.edu/ "PROBCONS").
