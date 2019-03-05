function [K,K_all] = gram_motion(K_path, D_path, sig, r_data)
% Computing the similarity based Kernel-Matrix using data distances

if ~isempty(K_path)
    load(K_path)
    if ~exist('K')
        error('Please save your kenel as variable "K"');
    end
    if nargout>1 && ~exist('K_all')
        error('Please also save your feature-wise kenel as variable "K_all"');
    end
    return
end

if isempty(D_path)
    error('Please provide the path to the saved distance matrix as "D_path"');
end

if ~exist('sig')
    sig=0;
end

%%
load(D_path);
if nargout>1 && ~exist('alldistances_feats')
    error('Please provide dimension-wise distance matrix as "alldistances_feats"');
end
if isempty(r_data)
    r_data=1:size(G,1);
end

G=alldistances(r_data,r_data);
G=G.^(0.5);     % apply this if distances are calculated by the DTW_calc.m
if nargout>1
    Gall=alldistances_feats(r_data,r_data,:);
    Gall=Gall.^(0.5);  % apply this if distances are calculated by the DTW_calc.m
end



if isempty(sig) || sig==0
    sig=mean(mean(G));
end
K = exp(-(G.^2 / (2 * sig.^2)));

if nargout>1
    sig=squeeze(mean(mean(Gall)));
    sig(:)=sig;
    for i=1:size(Gall,3)
        %                 param4(i)=squeeze(mean(mean(Gall(:,:,i))));
        K_all(:,:,i) = exp(-(Gall(:,:,i).^2 / (2 * sig(i)^2)));
    end
    
    if find(sig==0)
        Const_dims=find(sig==0)'
        fprintf('Constant dimensions are removed from the data!')
        pause(3);
    end    
end
