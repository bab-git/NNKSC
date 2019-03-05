function TSNE_visual(results,params)
% Approximate 2D visualization of:
% - data points 
% - the dictionary formation
% - sparse codings
% - supervised mapping (for LC-NNKSC)



show_tsne=0;
no_dims=2 ; max_iter=10000;
labels=results.list_labels;
r_tr=results.r_tr;
list_labels_tr=labels(r_tr);
if params.train_ratio<1
    r_ts=results.r_ts;
    list_labels_ts=labels(r_ts);
    r_val=results.r_val;
    list_labels_val=labels(r_val);
end

Adc=results.Adc;
Dic=dic_clean(Adc);
load (params.D_path,'alldistances')
% load (params.D_path,'alldistances','alldistances_feats')

D=alldistances;

n = size(D,1);
P = eye(n) - repmat(1/n,n,n);
S = P * (-.5 * D .* D) * P;
S_cor=eig_correct(S,'zero',0);
DisSq = bsxfun(@plus, bsxfun(@plus, -2.*S_cor, diag(S_cor)), diag(S_cor)');
D_new=sqrt(DisSq);
%% t-SNE
perplexity = floor(log(n)^2);
% perplexity =30;

MappedData = tsne_d(D_new, labels, no_dims, perplexity,max_iter,show_tsne);

plotData( MappedData, labels , 1);
title('APPROXIMATE 2-D visualization of original data');
%%  Prototypes and their constructions
prot_y=[];
for i_d=1:size(Dic,2)
    temp_i=find(Dic(:,i_d));
    prot_y{i_d}=[temp_i];
end

y_prot=[];
for i_n=1:size(Dic,1)
    temp_i=find(Dic(i_n,:));
    y_prot{i_n}=[temp_i];
end

% plotData_MP( MappedData,r_tr, labels, 1, prot_y);
plotData_MP2( MappedData,r_tr, labels, 1, y_prot);

%% TSNE sparse codes
if params.train_ratio==1
    X=dic_clean(results.X,0.1);
    labelsx=list_labels_tr;
else
    X=dic_clean(results.X,0.1);
    X_ts=dic_clean([results.X_ts],0.1);
    X_v=dic_clean([results.X_val],0.1);
    X=[X X_ts X_v];
    labelsx=[list_labels_tr list_labels_ts list_labels_val];
end

% for i_c=1:max(list_labels_tr)
%    i_cl=find(list_labels_tr==3)
% end

D=squareform(pdist(X'));
% labels=list_labels;

n = size(D,1);
P = eye(n) - repmat(1/n,n,n);
S = P * (-.5 * D .* D) * P;

% S=abs(X_ker)'*abs(X_ker);

% S=abs(fullX')+abs(fullX);
% S_corrected=eig_correct(S,params.kernel_modif,params.eig_correct_thresh);
% S_corrected=eig_correct(S,params.kernel_modif,0);
S_cor=eig_correct(S,'zero',0);
% clear S
% the way back (K -> D)
% S_cor=S;
DisSq = bsxfun(@plus, bsxfun(@plus, -2.*S_cor, diag(S_cor)), diag(S_cor)');
D_new=sqrt(DisSq);

N=size(D_new,2);
perplexity = floor(log(N)^2);
show_tsne=0;
D=D_new; no_dims=2 ; max_iter=10000; %perplexity=30;

MappedDataX  = tsne_d(D, labelsx, no_dims, perplexity,max_iter,show_tsne);

plotData( MappedDataX, labelsx , 1);
% plotData( MappedData(r_tr,:), list_labels_tr , 1);
titl={'2-D visualization of the sparse codings (X)';
    '**Sub-class Formations**'};
title(titl);

if strcmp(params.algorithm,'NNKSC')
    return
end
titl={'2-D visualization of the sparse codings (X)';
    '**Intermediate Mapping**'};
title(titl);
%% TSNE sparse codes  LC-kksc
if params.train_ratio==1
    X=dic_clean(results.X,0.1);
    labelsx=list_labels_tr;
    % r_ts=results.r_ts;
    % X=dic_clean([results.X_ts results.X_val],0.1);
    % labels=[list_labels_ts list_labels_val];
else
    X=dic_clean(results.X,0.1);
    X_ts=dic_clean([results.X_ts],0.1);
    X_v=dic_clean([results.X_val],0.1);
    X=[X X_ts X_v];
    labelsx=[list_labels_tr list_labels_ts list_labels_val];
end

Dic=dic_clean(Adc);
H_tr=0;
for i=1:length(r_tr)
    H_tr(list_labels_tr(i),i)=1;
end
X_l=H_tr*Dic*X;
% for i_c=1:max(list_labels_tr)
%    i_cl=find(list_labels_tr==3)
% end

D=squareform(pdist(X_l'));
% labels=list_labels;

n = size(D,1);
P = eye(n) - repmat(1/n,n,n);
S = P * (-.5 * D .* D) * P;

% S=abs(X_ker)'*abs(X_ker);

% S=abs(fullX')+abs(fullX);
% S_corrected=eig_correct(S,params.kernel_modif,params.eig_correct_thresh);
% S_corrected=eig_correct(S,params.kernel_modif,0);
S_cor=eig_correct(S,'zero',0);
% clear S
% the way back (K -> D)
% S_cor=S;
DisSq = bsxfun(@plus, bsxfun(@plus, -2.*S_cor, diag(S_cor)), diag(S_cor)');
D_new=sqrt(DisSq);

N=size(D_new,2);
perplexity = floor(log(N)^2);
show_tsne=0;
D=D_new; no_dims=2 ; max_iter=10000; %perplexity=30;

MappedDataX  = tsne_d(D, labelsx, no_dims, perplexity,max_iter,show_tsne);

plotData( MappedDataX, labelsx , 1);
% plotData( MappedData(r_tr,:), list_labels_tr , 1);
titl={'2-D visualization of the Supervised Mapping (H*A*X)';
    '**Sub-class Formations**'};
title(titl);
%% map protytopes 
n_D=size(Dic,2);
Kyy=results.Kyy;
Dx=Kyy(1)-Kyy;
D_xd=Kyy(1)-Kyy*Dic;
D_dd=Kyy(1)-Dic'*Kyy*Dic;
D=[[Dx D_xd];[D_xd' D_dd]];
labels=list_labels_tr;
ll=H_tr*Dic;
[v,l_d]=max(ll);
labels_MP_data=[labels l_d];

labels_all=list_labels_tr;


n = size(D,1);
P = eye(n) - repmat(1/n,n,n);
S = P * (-.5 * D .* D) * P;

S_cor=eig_correct(S,'zero',0);
DisSq = bsxfun(@plus, bsxfun(@plus, -2.*S_cor, diag(S_cor)), diag(S_cor)');
D_new=sqrt(DisSq);

N=size(D_new,2);
% perplexity = floor(log(N)^2);
perplexity =ceil(N/max(labels_all));
show_tsne=0;
no_dims=2 ; max_iter=10000; %perplexity=30; 
show_num=0;
markerSize=10;
font=14;

MappedData_prots_data  = tsne_d(D_new, labels_MP_data, no_dims, perplexity,max_iter,show_tsne);
plot_MP_data( MappedData_prots_data, labels_MP_data , n_D,1,show_num,font,markerSize,labels_all);
titl={'2-D visualization of the dictionary prototypes and the data'};
title(titl);
