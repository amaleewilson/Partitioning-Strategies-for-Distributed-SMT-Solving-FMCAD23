# Partitioning-Strategies-for-Distributed-SMT-Solving-FMCAD23
Data and guide for reproducing figures for the following paper https://arxiv.org/abs/2306.05854 

As mentioned in the paper, all of these experiments were run on a cluster with 
26 nodes running Ubuntu 20.04 LTS, each with 128 GB of RAM, and two Intel Xeon 
CPU E5-2620 v4 CPUs with 8 cores per CPU.

cvc5 and OpenSMT2 were used for partitioning, but only cvc5 was used to solve 
the partitions that were made by both solvers. 
cvc5 was allowed up to 60 seconds to partition a problem and up to 20 minutes to 
solve it.  OpenSMT2 had a partitioning time limit of up to 21 minutes, and it 
was also allowed up to 20 minutes of solving time. However, any problem that 
took more than 20 minutes in total (partitioning plus solving time) was counted 
as a timeout. 

Scrambling also had no timeout for scrambling time, but like partitioning, 
total time (scrambling time plus solving time) must not exceed 20 minutes. 
Otherwise, it is counted as a timeout. 

cvc5 version: 
https://github.com/amaleewilson/cvc5/tree/fmcad-partitioning-2
commit: 4b3bcf90cb7304896f0d3ae47836ba1d003ef780

OpenSMT2 version: 
https://github.com/usi-verification-and-security/opensmt
commit: 1cdc185a7e8ab0030d69259072bb75fdad7cf5c2

The SMT-COMP scrambler was used for scrambling:
https://github.com/SMT-COMP/scrambler
commit 2f2dbcd69d98894031c6359add0a898cd071bd98

data/partitioning_or_scrambling contains the time to partition or scramble the
problems, with file names corresponding to the strategy that was used. 

data/runs contains child folders for each strategy, with files named after the 
benchmark being solved.

data/logs contains any information realted to partitioning failures, which was
a rare occurrence, but sometimes a given strategy may work for up to 64 
partitions but then fail for 128 partitions, for example. 

benchmark_list.py contains the names of all of the benchmarks from SMT-COMP, 
QF_LRA, and QF_UF that were used. Note that the SMT-COMP benchmarks include
the QF_LIA, QF_RDL, and QF_IDL benchmarks. 

`python3 process_data.py`` will generate all figures from the paper in the 
figures directory. 