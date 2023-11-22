# Introduction

Finding the longest common subsequence (LCS) of multiple strings is an NP-hard problem, with many applications in the areas of bioinformatics and computational genomics. Although significant efforts have been made to address the problem and its special cases, the increasing complexity and size of biological data require more efficient methods applicable to an arbitrary number of strings. Here we present several algorithms for the general case of multiple LCS (or MLCS) problem, i.e., finding an LCS of any number of strings, and its parallel realization. These algorithms are based on the dominant point approach and employs a fast divide-and- conquer technique to compute the dominant points. 

# Quick-DP and Quick-DPPAR

Our algorithm Quick-DP finds exact LCS. When applied to a case of three strings, Quick-DP demonstrates the same performance as the fastest existing MLCS algorithm designed for that specific case. When applied to more than three strings, it is significantly faster than the best existing sequential methods, reaching up to 2-3 orders of magnitude faster speed on large-size problems. Quick-DPPAR is an efficient parallel implementation of Quick-DP. Evaluating the parallel algorithm on a benchmark set of both random and biological sequences reveals a near-linear speedup with respect to the sequential algorithm. The implementation of the both algoirthms are provided in mlcsparallel.cpp.

Citation
------------
[1] Q. Wang, D. Korkin, and Y. Shang (2011) A Fast Multiple Longest Common Subsequence (MLCS) Algorithm. IEEE Transactions on Knowledge and Data Engineering, 23(3):321-34.

[2] Q. Wang, D. Korkin, and Y. Shang (2009) Efficient Dominant Point Algorithms for the Multiple Longest Common Subsequence (MLCS) Problem. Intl. Joint Conf. on Artificial Intelligence (IJCAI).

 
# MLCS-Astar

With MLCS being an NP-hard problem, we also developed heuristic algorithms that are applicable to an arbitrarily large number of strings. MLCS-A* is a variant of the A* algorithm. It maximizes a new heuristic estimate of the LCS in each search step so that the longest common subsequence can be found. As a natural extension of MLCS-A*, a fast algorithm, MLCS-APP, was also developed to deal with large volume of biological data for which finding a LCS within reasonable time is impossible. The benchmark test shows that MLCS-APP is able to extract common subsequences close to the optimal ones and that MLCS-APP significantly outperforms existing heuristic approaches. When applied to 8 protein domain families, MLCS-APP produced more accurate results than existing multiple sequence alignment methods. The implementation of the algoirthm is provided in heuristic.cpp. For detailed description of this algorithm, please read the corresponding paper below.

Citation
------------
Q. Wang, M. Pan, Y. Shang, and D. Korkin (2010) A fast heuristic search algorithm for finding the longest common subsequence of multiple strings.  Twenty-Fourth AAAI Conference on Artificial Intelligence (AAAI-2010). 
