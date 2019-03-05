function [IP,D_intakes,D_intakes_ind]=D_interpretability(results,thresh)
% The code checks how interpretabile is the learned dictionary.

% OUTPUTS:
% IP = Interpretability measure for dictionary elements. 100% --> the
     % dictionary element uses only one class of data
% D_intakes_ind: A matrix which shows list of the classes that each
     % dictionary element uses
% D_intakes = More detailed information on D_intakes and reelvant weights
     % in the dictionary elements.

%INPUTS:

%results= results of the NNKSC or LC-NNKSC algorithms.
%thresh= a thresh based on which dictionary elements can be pruned.


Adc=results.Adc;
r_tr=results.r_tr;
list_labels_tr=results.list_labels(r_tr);
H_tr=results.H_tr;

if ~exist('thresh')
    thresh=0.25;
end
MP=size(Adc,2);
%%
Adc_cln=dic_clean(Adc,thresh);

D_intakes=[];
D_intakes_ind=zeros(MP,size(H_tr,1));

for i_d=1:size(Adc,2)
    D_intakes{i_d}=[list_labels_tr(find(abs(Adc_cln(:,i_d))));Adc_cln(find(abs(Adc_cln(:,i_d))),i_d)'];
    temp=unique(D_intakes{i_d}(1,:));
    D_intakes_ind(i_d,find(temp))=temp;
end
if isempty(find(D_intakes_ind(:,2:end)))
    disp('***** No interconnction of different classes in the dictionary *****')
else
    i_ov=find(D_intakes_ind(:,2))';
    fprintf('Dictionary element %d has class-overlappings\n', i_ov);
end
%%

atom_share_total=zeros(size(H_tr,1));
for i_d=1:size(Adc_cln,2)    
    temp_D=abs(Adc_cln(:,i_d));
%     temp_D(find((temp_D/max(temp_D)*100)<thresh))=0;
    atom_share=H_tr*abs(temp_D);
    IP(i_d)=max(atom_share)/sum(atom_share)*100;
    class_share(:,i_d)=atom_share/max(atom_share);
    [val, i_class]=max(atom_share);
%     atom_share=atom_share/max(atom_share);
%     atom_share_total(i_class,:)=atom_share_total(i_class,:)+atom_share';    
end
% atom_share_total=atom_share_total./repmat(max(atom_share_total,[],2),1,size(H_tr,1))*100;
temp=sort(atom_share_total,2,'descend');

[val bad_i]=min(IP);
sbad=val;
[val good_i]=max(IP);
sgood=val;
% IP_all=diag(atom_share_total);
%%
figure
stem(class_share(:,bad_i))
xlabel('classes')
ylabel('Intake Precentage (%)')
grid on;
title(sprintf('Worst interpretable dictionary element = %d , IP= %3.2f%%',bad_i,sbad))
figure
stem(class_share(:,good_i))
xlabel('classes')
ylabel('Intake Precentage (%)')
grid on;
title(sprintf('Best interpretable dictionary element = %d , IP= %3.2f%%',good_i,sgood))

fprintf('Best interpretability for dictionary elements (bIP)= %3.2f%% \n', sgood)
fprintf('Worst interpretability for dictionary elements (wIP)= %3.2f%% \n', sbad)
fprintf('Average interpretability for dictionary elements (aIP)= %3.2f%% \n', mean(IP))