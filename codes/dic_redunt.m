function [Ad_red,rec_err,X]=dic_redunt(Ad,Kyy,params,X0,rate,thresh)
% Cleaning the redundancies inside the dictionary

if strcmp(params.algorithm,'LC-NNKSC')   % LC-NNKSC mapping is more sensetive to the dictionary structure
    rate=1;
end
    
Ad0=Ad;
i_z0=find(sum(Ad,1)==0);

i_zx=find(sum(X0,2)==0)';
i_z=[i_z0 i_zx];
Ad(:,i_zx)=[];
Ad0_red=Ad;

if ~exist('thresh')
    thresh=0.25;
end
Ad_cln=dic_clean(Ad,thresh);
% Ad(Ad_cln==0)=0;

if ~exist('rate')
    rate=0.95;
end

%===== reduction
rec_err0=PMSE_kern(Ad0,X0,Kyy,Kyy,Kyy);
% i_ch=0;
for rate2=rate:0.005:1
    [Ad_red,X,rec_err]=dic_redunt_sub(Ad,Kyy,rate2,params);
    err_change=(rec_err-rec_err0)/rec_err0*100;
    if err_change<=10
        fprintf('prunning error= %3.2f%% , similarity threshold= %3.2f%% \n', err_change,rate2*100)
        return    
    end
end


% ========
function [Ad,X,rec_err]=dic_redunt_sub(Ad,Kyy,rate,params)
T=params.T;
SR=params.self_rec;
MP_n=size(Ad,2);
id=1;
while id<= MP_n
    cand=sparse(Ad(:,id));
    auxD=Ad;
    auxD(:,id)=0;
    simi = abs(cand'*Kyy*auxD) ;
    i_f=find(simi>rate);
    if ~isempty(i_f)
        [v ix]=max(simi(i_f));
        i_max=i_f(ix);
        Ad(:,id)=[];
    else
        id=id+1;
    end
    MP_n=size(Ad,2);
end
%==calculate X and error
X=nnkomp_all(Ad,Kyy,Kyy,Kyy,T,SR);
rec_err=PMSE_kern(Ad,X,Kyy,Kyy,Kyy);