% Usage: [mynorm] = Knorms(A,K)
% compute norm of columns in feature space
% given coefficients matrix A and Gram matrix
function [mynorm] = Knorms(A,K)

mynorm = zeros(size(A,2),1);

for i=1:size(A,2) 
    mynorm(i)  = sqrt(max(A(:,i)'*K*A(:,i),0)) ; 
end

end