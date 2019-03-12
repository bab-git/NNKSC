# None-Negative Kernel Sparse Coding (NNKSC)
NNKSC is a kernel-based sparse coding and dictionary learning algorithm which enforces non-negativity constraints on the dictionary and the sparse codes. As a result, the learned dictionary atoms and the sparse encodings are more interpretable regarding the semantic characteristics.
## Using NNKSC
- The NNKSC_main.m contain a test run of the NNKSC algorithm on 3 databases.
- In order to use NNKSC in the supervised setting, (LC-NNKSC) the parameter LC_betta should be tuned to a non-zero weight.
- More instruction is provided in the read_me.txt file.


## The paper
Hosseini, Babak, et al. "Non-negative kernel sparse coding for the analysis of motion data.", *International Conference on Artificial Neural Networks (ICANN)*, Springer, Cham, 2016.
