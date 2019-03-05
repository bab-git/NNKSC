# None-Negative Kernel Sparse Coding (NNKSC)
NNKSC is a kernel-based sparse coding and dictionary learning algorithm which enforces non-negativity constraints on the dictonary and the sparse codes. As a result, the learned dictionary atoms and the sparse encodings are more interpretable regarding the semantic charactristics.
## Using NNKSC
The NNKSC_main.m contain a test run of the NNKSC algorithm on 3 databases.
In order to use NNKSC in the supervised setting, (LC-NNKSC) the parameter LC_betta should be tuned to a non-zero weight.
