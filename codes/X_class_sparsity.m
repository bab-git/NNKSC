function [sparse_cls]=X_class_sparsity(results,thresh)
labels=results.list_labels_tr;
X=full(results.X);
if ~exist('thresh')
    thresh=0.05;
end
sbig=0;
ssmall=inf;
for i=1:max(labels)
    class_i=find(labels==i);
    sparse_perf=sum(abs(X(:,class_i)),2);
    ss_norm=sparse_perf/max(sparse_perf);
    sparse_perf(ss_norm<thresh)=0;
    num_sparse=numel(find(sparse_perf));
    sum_temp=mean(ss_norm);
    ss2=(ss(1)-ss(2))/ss(1)*100;
    if num_sparse > sbig
        i_big=i;
        sparse_bad=sparse_perf;
        sbig=num_sparse;
    elseif num_sparse < ssmall
        i_small=i;
        sparse_good=sparse_perf;
        ssmall=num_sparse;
    end
    sparse_cls(i)=sum(sparse_perf>0);
end

figure
stem(sparse_good)
grid on;
title(strcat('Best sparseness of X : ',num2str(ssmall),' used atoms  -  class Nr.',num2str(i_small)))

figure
stem(sparse_bad)
grid on;
title(strcat('Worst sparseness of X : ',num2str(sbig),' used atoms  - class Nr.',num2str(i_big)))

printf('Best class-wise sparsity (bSP)= %d \n', sum(sparse_good>0))
printf('Worst class-wise sparsity (wSP)= %d \n', sum(sparse_bad>0))
printf('Average class-wise sparsity (aSP)= %3.2f \n', mean(sparse_cls))