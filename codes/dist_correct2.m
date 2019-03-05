function alldistances_new=dist_correct2(alldistances)
sx=max(max(max(alldistances)));
alldistances_new=zeros(size(alldistances));
for i=1:size(alldistances,3)
    temp=triu(alldistances(:,:,i))+tril(alldistances(:,:,i)');
    temp=temp-diag(diag(temp));
    alldistances_new(:,:,i)=temp;
end