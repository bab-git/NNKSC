function [MSE,MSE_abs] = PMSE_kern(A,X,Kyy,k_zz,k_zy,ix)
%  Computing the precentage MSE error for recosntruction of data
if ~exist('k_zz')  && ~exist('k_zy')
    k_zz=Kyy;
    k_zy=Kyy;    
end
if isempty(k_zz)  && isempty(k_zy)
    k_zz=Kyy;
    k_zy=Kyy;  
end
if exist('ix')
    Kyy=Kyy;
    k_zz=k_zz(ix,ix);
    k_zy=k_zy(ix,:);
    X=X(:,ix);
end

MSE_abs=trace(abs(-2*k_zy*A*X+k_zz+X'*A'*Kyy*A*X));
MSE=MSE_abs / trace(k_zz)*100;