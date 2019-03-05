function [class_accu, conf_mat, i_missed]=LC_classification(X,A,Htr,labels)
% U=Q*D;
W=Htr*A;
% X_class=X(:,1);
% class_num=max(labels);
labels_test=W*X;
class_res=abs(ones(size(Htr,1),size(X,2))-labels_test);
[val i_class]=min(class_res,[],1);
% (i_class==labels)
class_accu=sum(i_class==labels)/length(labels)*100;
conf_mat=zeros(max(labels));
i_missed=[];
for i_y=1:size(X,2)
    conf_mat(labels(i_y),i_class(i_y))=conf_mat(labels(i_y),i_class(i_y))+1;
    if i_class(i_y)~=labels(i_y)
        i_missed=[i_missed i_y];
    end
end