function [Kyy, P, K_modif_perc]=eig_correct(Kyy,ker_modif,rate)
P=0;
if isempty(rate) || rate==0
    rate=0.5;
end
if abs(min(eig(Kyy))/max(eig(Kyy))*100)>rate && min(eig(Kyy))<0
    K0=Kyy;
    [V,eigK]=eig(K0);
    eigK2=(sqrt(eigK.^2));
    K2=V*eigK2*V^-1;
    eigK3=eigK;
    eigK3(eigK3<0)=0;
    K3=V*eigK3*V^-1;
    switch ker_modif
        case 'zero'
            Kyy=K3;  % make negative eigs to zero
            M=eigK3~=0;
            P=V*M*V';
        case 'abs'
            Kyy=K2;  % make negative eigs to abs(eig)
            M=sign(eigK);
            P=V*M*V';
    end
    K_modif_perc=norm(K0-Kyy)/norm(K0)*100;
    disp('--Kernel is modified regarding its negative eigen-values')
end