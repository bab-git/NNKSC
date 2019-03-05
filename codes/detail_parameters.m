%% detailed parameters --> not have direct effect on the performance
params.self_rec=1;   % 1: lets each xi use dictionary elements aj in which also xi has a part
                     % 0: no self-reconstruction ==> leads to more sub-class devisions
params.red_thresh=.25;   % the threshold chosen to eliminate entires of dictionary
params.sim_thresh=0.99;  % the threshold chosen to eliminate similar dictionary elements
params.i_term=10;     % number of iterations to decide for stoping of algorithm
params.loop_a=1;     % Number of internal iterations for updating each ai column
params.eig_correct_type='zero';   % how to correct negative eigen-values in the kernel --> 'zero'=cliping , 'abs'=flipping
params.eig_correct_thresh=0.1;