
This code is related to the publication

Hosseini B, Hülsmann F, Botsch M, Hammer B, 
"Non-Negative Kernel Sparse Coding for the Analysis of Motion Data" 
Artificial Neural Networks and Machine Learning – ICANN 2016


================ General description
The algorithms optimizes the following spares coding framework:
<A,X>=argmin ||Phi(Y)-Phi(Y)*A*X||
s.t.  |x|_0 <T , |a|_0 <T 
      x >=0 , a >=0
It uses the similarity kernel of motion data points based on Dynamic Time Warping.
It learns a dictionary based representation for the data showing 
how data points can be reconstructed using the learned dictionary prototypes.

In the supervised setting (LC-NNKSC), the annotations information is used to make the 
represetation for meaningful regarding the supervised information.

NNKSC_main.m:      The code for "Non-Negarive Kernel Sparse Coding" algorithm   (Unsupervised algorithm)
                   THe representation can reviel underlying hidden subcategories without using the labeling information

LC_NNKSC_main.m:   The code for "LC-Non-Negarive Kernel Sparse Coding" algorithm   (Supervised algorithm)
                   it provides a supervised dictionary learning and mapping using the labeling information


** In case of having any problem or question regarding the code please contact:
bbkhosseini@gmail.com or bhosseini@techfak.uni-bielefeld.de

% Babak Hosseini
% CITEC, University Bielefeld 