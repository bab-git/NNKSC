function [X, res_x]= nnKomp(A, Kyy, Kzy,Kzz,T0,G,B)
% None-Negative Kernel Orthogonal Matching Pursuit
% solving the problem    min_x ||Phi(z)-Phi(Y)*A*x||_f    s.t  |x|_0<=T0
% Inputs:
%       A: dictionary matrix
%       Kyy: Gram-matrix of training data
%       Kzy: Gram-matrix of test-data respect to training data
%       Kzz: Gram-matrix of test-data
%       G=A'*Kyy*A;
%       B=A'*Kyy';
%       T0=A'*Kyy';
% Output:
%       X: sparse encoding for reconstruction of z
%       res_x: reconstruction error

% Hosseini and Hammer,"Non-Negative Kernel Sparse Coding for the Analysis of Motion Data", ICANN2016
%% ------------------------------
e_tol=1;
[~,N] = size(A);
X = sparse(zeros(N,1));
S = []; % positions indexes of components of s
R = 1:size(A,2); % positions indexes of components of s
R(sum(A,1)==0)=[];
res_phi = Kzy*A; % first r*d

x_est_pre=0;
res_x=100;
x2=zeros(N,1);r_1=0;
t=1;
while (t<=T0)
    [~,j]=min((res_phi(R)-1).^2);
    j=R(j);

    S = [S j];
    R=R(R~=j);
    S=sort(S);   
    Ai=A(:,S);    
    x_est = Klsqnonneg(Ai,Kyy,Kzy);    
    if exist('G')
        x2=X;x2(S)=x_est;        
        Gx2=G*x2;
        res_x2=(res_x/100*Kzz+r_1-x2'*Gx2)/Kzz*100;        
    else
        res_x2=PMSE_kern(Ai,x_est,Kyy,Kzz,Kzy);
    end
    
    
    if (res_x2-res_x) < -res_x*0.001
        res_x=res_x2;
        t=t+1;
    else
        S=S(S~=j);
        Ai=A(:,S);
        if isempty(Ai)
            Ai=zeros(size(Ai,1),1);
        end
        x_est=x_est_pre;
        if exist('G')
            x2=X;x2(S)=x_est;
            Gx2=Gx2_0;
        end
        
    end
    
    if exist('G')
        res_phi=B'-Gx2';  %r'*D
        r_1=x2'*G*x2;
        Gx2_0=Gx2;
    else
        vs = Ai*x_est;
        res_phi=(Kzy-(vs)'*Kyy)*A; %r'*D
    end
    
    x_est_pre=x_est;
    
    if norm(res_x) < e_tol  || isempty(R) 
        break
    end
    
    
end
X(S)=x_est;
