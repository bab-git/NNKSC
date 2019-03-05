function [ai, fobj]= my_nnsign_A(A, xj, K, XtX, xd, Ej,ro, params,i_a, nn_flag, kmp,idy)
% L, Hl, H, params.lambda,params.lambda)
if isfield(params,'cl_T_a')
    T0=params.cl_T_a;
    if T0==0
        T0=size(A,1);
    end
elseif exist('kmp')
    T0=kmp;
else
    T0=size(A,1);    
end

if ~exist('idy')
    idy=[];
end

if ~exist('nn_flag')
    nn_flag=0;
end
% K=Kyy;
% A=D_kksvd;
%


% if exist('no_lam')
if T0<size(A,1)
    lambda=0;
else
    lambda=params.lambda;
end
% end

if isempty(ro)
    bet=0;
else
    bet=params.cl_betha;
end
gam=params.a_gamma;

% lambda=0.027
% gam=1


% qp_solve=params.cl_qp_solve;
%==== constructing H and c
N=size(A,1);
% a0=A(:,i_a);
% X0=X_ker0;
% xj=X0(:,i_a);


Hi=ones(N,1)*xd'-XtX;
% Hi=rand(size(Hi));
Hi=(Hi'+Hi)/2;

% H_all=H;
if bet==0
    H_all=(xj*xj'*K)+gam*Hi;
else
    H_all=(xj*xj'*K)+bet*(K*ro(i_a,i_a))+gam*Hi;
end
% H_all=alp*H +bet*(H*L(i_a,i_a))+gam*(Hl(i_a,i_a));
% H_all=H +(H*L(i_a,i_a))+(Hl(i_a,i_a));
% H_all=0.5*(H_all+H_all');
H_all=2*H_all;
% eig(H)
eh=eig(H_all);
if abs(min(eh)/max(eh))*100>.1 && min(eh)<0
    negative_eig_ratio=abs(min(eh)/max(eh))*100;
    %     error('min(eh)<0')
    %     display('min(eh)<0')
%     qp_solve=1;
end
A_i=full(A);
A_i(:,i_a)=0;
if bet==0
    c_all=(-2*xj*Ej'*K);
else
    c_all=(-2*xj*Ej'*K) + bet*(K*A_i*ro(:,i_a))';
end
% res0=0.5*x0'*H_all*x0+c_all*x0+K(i_a,i_a);
%% solving feature - sign
% [ai,fobj]= NNl1ls_featuresign(H_all, c_all', lambda, 0);
params.cl_FSS=1;
if params.cl_FSS
    if max(abs(c_all))<lambda
        lambda=max(abs(c_all));
    end
    
        [ai,fobj]= NNl1ls_featuresign2(H_all, c_all', 0*lambda, 0,T0,idy);
    
else
    [ai,fobj]= NNl1ls_QP(H_all, c_all', lambda);
end

while sum(ai)==0 && ~params.cl_FSS
    lambda=lambda*0.75;
    if params.cl_FSS
        [ai,fobj]= NNl1ls_featuresign(H_all, c_all', lambda, 1,qp_solve);
    else
        [ai,fobj]= NNl1ls_QP(H_all, c_all', lambda);
    end
    if lambda<0.001
        %         display('zero ai even with small lambda !!!')
        break
    end
    %     printf('reduced lambda_a= %3.3f',lambda)
end

% fobj_x(x_new)
% fobj_x(ai)
% list_labels_tr(find(x_new))
% list_labels_tr(find(ai))
% numel(find(x_new))
% numel(find(ai))


% ai=ai/sqrt(ai'*K*ai);
obj=@(ai) 0.5*ai'*H_all*ai+c_all*ai;
fobj=obj(ai);
obj_a=ai'*Hi*ai;
return
%%
% i=randperm(58,1)
x1=X_ker1(:,i);x=X_ker(:,i);x0=X_ker0(:,i);
res0=0.5*x0'*H_all*x0+c_all*x0+K(i_a,i_a)
res1=0.5*x1'*H_all*x1+c_all*x1+K(i_a,i_a)
res=0.5*x'*H_all*x+c_all*x+K(i_a,i_a)
% xf=full(x);
% komp_res=((-2*K(i_a,:)*A*xf+K(i_a,i_a)+xf'*H*xf))
% X2=X_i;X2(:,i_a)=x;
% trace(X2'*H*X2*L)
% trace(X2'*X2*Hl)