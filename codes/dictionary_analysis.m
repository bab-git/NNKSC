%% ==== clear uneeded files + ordering the D
dic_red_thresh=params.red_thresh;
Adc_cleaned=Adc./repmat(max(abs(Adc)),[size(Adc,1) 1]);
Adc_cleaned(abs(Adc_cleaned)<dic_red_thresh)=0;
Adc_cln=Adc; Adc_cln(Adc_cleaned==0)=0;
% X_ker_ts=[];X_ker=[];X_ker_val=[];
% param2=[];
% for i_y=1:numel(r_tr)
%     %         save KOMP_state
%     Adc_temp=Adc;
%     id=find(Adc(i_y,:));
%     id2=sum(Adc(:,id)>0);
%     if find(id2==1)
%         id2;
%     end
%     Adc_temp(:,id(find(id2==1)))=0;
%     [X_ker(:,i_y), res_x]= my_nnKomp_1(Adc_temp, Kyy, Kyy(i_y,:),Kyy(i_y,i_y),T0);
% end
% if params.train_ratio<1
%     for i_y=1:numel(r_ts)
%         if params.OMP_res_1==0
%             [X_ker_ts(:,i_y), res_x]= my_nnKomp(Adc_nrm, K0, Kzy(i_y,:),Kzz(i_y,i_y),T0);
%         else
%             [X_ker_ts(:,i_y), res_x]= my_nnKomp_1(Adc_nrm, K0, Kzy(i_y,:),Kzz(i_y,i_y),T0);
%         end
%     end
%     fullX_ts=full(X_ker_ts);
%     for i_y=1:numel(r_val)
%         if params.OMP_res_1==0
%             [X_ker_val(:,i_y), res_x]= my_nnKomp(Adc_nrm, K0, Khy(i_y,:),Khh(i_y,i_y),T0);
%         else
%             [X_ker_val(:,i_y), res_x]= my_nnKomp_1(Adc_nrm, K0, Khy(i_y,:),Khh(i_y,i_y),T0);
%         end
%     end
%     fullX_val=full(X_ker_val);
%
% end
%
% MSEP_tr=error_kern(Adc,X_ker,Kyy,Kyy,Kyy)
% if params.test_ksvd==1
%     MSEP_test=error_kern(Adc,X_ker_ts,Kyy,Kzz,Kzy)
%     MSEP_val=error_kern(Adc,X_ker_val,Kyy,Khh,Khy)
% end
% fullX=full(X_ker);
% % A = Adc ;
% % MSEP=mean(komp_res)*100
%
% % end
%% similar classes check (vario)
for i=1:size(Kyy,1)
    i_dom=[];
    i_minor=[];
    %             fi{i}=[fi{i} fullX(:,find(i_list_var==i2))];
    %     fj=fullX(:,find(i_list_var==j));
    %     fk=fullX(:,find(i_list_var==k));
    %             for i3=1:size(fi{i},2)
    temp=fullX(:,i);
    temp=temp/max(temp)*100;
    i_dom=[i_dom find(temp>30)'];
    i_minor=[i_minor find((temp<30).*(temp>1))'];
    %             end
    
    list_dom{i}=unique(i_dom);
    %         list_dom.(Y_segments{i})=unique(i_dom);
    list_minor{i}=unique(i_minor);
end
%% dictionary sparsity
thresh=10; % Threshold for dominancy of signals in D
D_intakes=[];
D_intakes_ind=zeros(MP,max(list_labels));

for i_d=1:size(Adc,2)
    D_intakes{i_d}=[list_labels_tr(find(abs(Adc_cleaned(:,i_d))>thresh));Adc_cleaned(find(abs(Adc_cleaned(:,i_d))>thresh),i_d)'];
    temp=unique(D_intakes{i_d}(1,:));
    D_intakes_ind(i_d,find(temp))=temp;
end
if isempty(find(D_intakes_ind(:,2:end)))
    display('***** No interconnction of different classes in the dictionary *****')
end
% break