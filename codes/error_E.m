function norm_Ky=error_E(E,Ky,A,x)
% norm_Ky=0;norm_E=0;
% for i_n=1:size(E,2)
%     norm_Ky=norm_Ky+E(:,i_n)'*Kyy*E(:,i_n)-(2*E(:,i_n)'-x(i_n)'*A')*KyA*x(i_n);
%     norm_E=norm_E+E(:,i_n)'*Kyy*E(:,i_n);
% end
% norm_Ky=norm_Ky/norm_E*100;

norm_Ky=trace((E-A*x)'*Ky*(E-A*x))/trace(E'*Ky*E)*100;