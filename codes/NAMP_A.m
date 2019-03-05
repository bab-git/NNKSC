function [ai, fobj]= NAMP_A(A, xj, K, XtX, xd, Ej, params,idy)
% ai0=sparse(A(:,i_a));
if isfield(params,'T_a')
    T0=params.T_a;
    if T0==0
        T0=size(Ej,1);
    end
else
    T0=size(Ej,1);
end

if ~exist('idy')
    idy=[];
end

gam=params.a_gamma;

%==== constructing H and c
N=size(A,1);

Hi=ones(N,1)*xd'-XtX;
Hi=(Hi'+Hi)/2;

H_all=(xj*xj'*K)+gam*Hi;

H_all=2*H_all;

c_all=(-2*xj*Ej'*K);
%% solving feature - sign
[ai,fobj]= NAMP(H_all, c_all', 0, 0,T0,idy);
obj=@(ai) 0.5*ai'*H_all*ai+c_all*ai;
fobj=obj(ai);
obj_a=ai'*Hi*ai;