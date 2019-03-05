function results=NNKSC_func(params)
% The algorithms optimizes the following spares coding framework:
% <A,X>=argmin ||Phi(Y)-Phi(Y)*A*X||
% s.t.  |x|_0 <T , |a|_0 <T
%       x >=0 , a >=0
% It uses the similarity kernel of motion data points based on Dynamic Time Warping.
% It learns a dictionary based representation for the data showing
% how data points can be reconstructed using the learned dictionary prototypes.
%
% In the supervised setting (LC-NNKSC), the annotations information is used to make the
% represetation for meaningful regarding the supervised information.

list_labels=params.list_labels;
T0 = params.T;
iternum = params.iter_end ;
MP = params.MP;
K_path=params.K_path;
D_path=params.D_path;
i_term=params.i_term;
eig_correct_type=params.eig_correct_type;
eig_correct_thresh=params.eig_correct_thresh;

if strcmp(params.algorithm,'LC-NNKSC')
    lcNN=1;
else
    lcNN=0;
end

% if lcNN==1 && params.train_ratio==1
%     params.train_ratio=0.75;
% end

n_data=length(list_labels);
if params.train_ratio<1    
    [r_tr,r_rest]=makesplits_CV(list_labels,params.train_ratio,1,1,1);
    
    [r_ts,r_val]=makesplits_CV(list_labels(r_rest),0.5,1,1,1);
    r_ts=r_rest(r_ts);
    r_val=r_rest(r_val);
    
    
else
    r_tr=randperm(n_data);
    r_ts=[];r_val=[];
end
list_labels_tr=list_labels(r_tr);
list_labels_ts=list_labels(r_ts);
list_labels_val=list_labels(r_val);
n_tr=numel(r_tr);
%% initialize Dictionary
if params.self_rec==0
    params.dic_init=1;
end
Adc_init = zeros(n_tr, MP) ;
if MP>n_tr
    randid = randi(n_tr,ceil(MP/n_tr),MP);
    for i=1:MP
        Adc_init(randid(:,i),i) = 1;  % randomly initilize dictionary using input signals
    end
else
    randid = randperm(n_tr);
    for i=1:MP
        Adc_init(randid(i),i) = 1;  % randomly initilize dictionary using input signals
    end
end
% save temp_init
%% Computing the labeling kernel
H_tr=0;
for i=1:numel(r_tr)
    H_tr(list_labels_tr(i),i)=1;
end
%% Compute the kernel
Kyall= gram_motion(K_path, D_path, [],[r_tr r_ts r_val]); % compute Gram matrix
Kyall=eig_correct(Kyall,eig_correct_type,eig_correct_thresh);
mynormall = sqrt(diag(Kyall));
Kyall = Kyall./(mynormall*mynormall');

Kyy=Kyall(1:length(r_tr),1:length(r_tr));
Kyy0=Kyy;
if params.train_ratio<1
    Kzz=Kyall(length(r_tr)+[1:1:length(r_ts)],length(r_tr)+[1:1:length(r_ts)]);
    Kzy=Kyall(length(r_tr)+[1:1:length(r_ts)],1:length(r_tr));
    Khh=Kyall(length(r_tr)+length(r_ts)+[1:1:length(r_val)],length(r_tr)+length(r_ts)+[1:1:length(r_val)]);
    Khy=Kyall(length(r_tr)+length(r_ts)+[1:1:length(r_val)],1:length(r_tr));
end

i_tr=1:length(r_tr);
if lcNN
    Kh=H_tr'*H_tr;
    Kyy=Kyy0+params.LC_betta*Kh(i_tr,i_tr);%+LC_alfa*Kq(tr_r,tr_r);
end


% normalizing dictionary elements
for i=1:size(Adc_init,2)
    Ad=Adc_init(:,i);
    Adc_init(:,i)=Ad/sqrt(Ad'*Kyy*Ad);
end
Adc=Adc_init;

total_err = zeros(iternum,1) ;
XBIAS=0;
if params.display==1
    if lcNN
        fig1=figure;
        fig3=figure('Units','centimeters ','Position',[19-XBIAS 5 16 20],'PaperPositionMode','auto');
    else
        fig2=figure;
    end
    fig5=figure;
    %     fig3=figure('Units','centimeters ','Position',[5-XBIAS 3 13 13],'PaperPositionMode','auto');
    %     fig1=figure('Units','centimeters ','Position',[19-XBIAS 5 16 20],'PaperPositionMode','auto');
    %     fig5=figure('Units','centimeters ','Position',[35-XBIAS 3 13 13],'PaperPositionMode','auto');
end

iter=1;
X_total=[];
Adc_total=[];

err_recon=0;
err_recon_ts=0;
err_recon_val=0;
first_test=1;
const_perf=0;

%% ===== loop
while (iter<iternum)
    fprintf('%d...',iter);  %if (mod(iter,30)==0), fprintf('\n'); end
    
    X=nnkomp_all(Adc,Kyy,Kyy,Kyy,T0,params.self_rec);
    %     norm(find(X),'fro')
    
    err_obj(iter)=PMSE_kern(Adc,X,Kyy,Kyy,Kyy);
    err_recon(iter)=PMSE_kern(Adc,X,Kyy0,Kyy0,Kyy0);
    mapp_err(iter)= norm(H_tr-H_tr*Adc*X,'fro')^2/norm(H_tr,'fro')^2*100;
    if lcNN
        accur_tr_LC(iter)=LC_classification(X,Adc,H_tr,list_labels_tr);
    end
    %% test phase
    if params.train_ratio<1
        
        X_ts=nnkomp_all(Adc,Kyy0,Kzy,Kzz,T0);
        
        X_val=nnkomp_all(Adc,Kyy0,Khy,Khh,T0);
        
        err_recon_ts(iter)=PMSE_kern(Adc,X_ts,Kyy,Kzz,Kzy);
        err_recon_val(iter)=PMSE_kern(Adc,X_val,Kyy,Khh,Khy);
    end
    if lcNN && params.train_ratio<1
        accur_ts_LC(iter)=LC_classification(X_ts,Adc,H_tr,list_labels_ts);
        accur_val_LC(iter)=LC_classification(X_val,Adc,H_tr,list_labels_val);
    end
    
    X_total{iter}=X;
    Adc_total{iter}=Adc;
    A = eye(n_tr) - Adc*X ;
    proj_err = Knorms(A,Kyy) ; % norm of columns in feature space
    [~, LeastRepId] = sort(proj_err, 'descend') ;
    
    %% ======== figures
    if params.display==1 && iter>1
        if lcNN
            change_current_figure(fig3)
            subplot(3,1,1)
            hold off;
            plot([2:iter],err_obj(2:end),'b','marker','o'); title(sprintf('Whole objective error= %3.2f%%',err_obj(end)))
            grid on;
            
            subplot(3,1,2)
            hold off;
            plot([2:iter],err_recon(2:end),'b','marker','o'); title(sprintf('Reconstruction Obj error= %3.2f%%',err_recon(end)))
            grid on;
            hold on
            if params.train_ratio<1
                plot([2:iter],err_recon_ts(2:end),'r','marker','o'); %title('Reconstruction Error = MSE %')
                plot([2:iter],err_recon_val(2:end),'g','marker','o'); title(sprintf('Reconstruction Obj error= %3.2f%%',err_recon_val(end)))
                legend('train','test','validation')
            end
            
            subplot(3,1,3)
            hold off;
            plot([2:iter],mapp_err(2:end),'b','marker','o'); title(sprintf('Supervised Obj. error = %3.2f%%',mapp_err(end)))
            grid on;
            
            
            change_current_figure(fig1)
            hold off;
            plot([2:iter],accur_tr_LC(2:iter),'b','marker','o'); title(sprintf('LC-classification accuracy = %3.3f %',accur_tr_LC(end)))
            if params.train_ratio<1
                hold on;
                plot([2:iter],accur_ts_LC(2:iter),'r','marker','o');
                plot([2:iter],accur_val_LC(2:iter),'g','marker','o'); title(sprintf('LC-classification accuracy = %3.3f %',accur_val_LC(end)))
                
            end
            grid on;
        else
            change_current_figure(fig2)
            hold off;
            plot([2:iter],err_recon(2:end),'b','marker','o'); title(sprintf('Reconstruction error= %3.2f%%',err_recon(end)))
            grid on;
            hold on
            if params.train_ratio<1
                plot([2:iter],err_recon_ts(2:end),'r','marker','o'); %title('Reconstruction Error = MSE %')
                plot([2:iter],err_recon_val(2:end),'g','marker','o'); title(sprintf('Reconstruction error= %3.2f%%',err_recon_val(end)))
                legend('train','test','validation')
            end
        end
    end
    
    drawnow
    %% --- class-wise formation of dictionary
    thresh=params.red_thresh;
    Adc_cleaned=dic_clean(Adc,thresh);
    
    
    D_intakes=[];
    D_intakes_ind=zeros(MP,max(list_labels));
    for i_d=1:size(Adc,2)
        D_intakes{i_d}=[list_labels_tr(find(abs(Adc_cleaned(:,i_d))));Adc_cleaned(find(abs(Adc_cleaned(:,i_d))),i_d)'];
        temp=unique(D_intakes{i_d}(1,:));
        D_intakes_ind(i_d,find(temp))=temp;
    end
    shared_MP(iter)=numel(find(D_intakes_ind(:,2)));
    
    if params.display==1
        change_current_figure(fig5)
        hold on;grid on
        plot([1: iter],shared_MP,'marker','o'); title(sprintf('# Dictionary atoms with class-overlappings = %d',shared_MP(end)))
    end
    
    drawnow
    
    fprintf('..PMSE_tr= %3.2f%%',err_recon(iter));
    if lcNN
        fprintf('...Class_tr = %3.2f%%',accur_tr_LC(iter));
    end
    if (params.train_ratio<1)==1 % && params.display
        fprintf('..PMSE_ts = %3.2f%%',err_recon_ts(iter));
        fprintf('..PMSE val = %3.2f%%',err_recon_val(iter));
        if lcNN
            fprintf('...Class_ts = %3.2f%%',accur_ts_LC(iter));
            fprintf('...Class_val = %3.2f%%',accur_val_LC(iter));
        end
    end
    fprintf('\n');
    
    
    %% === terminaiton CHECK
    i_term=params.i_term;
    if params.train_ratio<1
        if lcNN
            performance_curve(iter)=100-accur_ts_LC(iter);
        else
            performance_curve(iter)=err_recon_ts(iter);
        end
    else
        if lcNN
            performance_curve(iter)=100-accur_tr_LC(iter);
        else
            performance_curve(iter)=err_recon(iter);
        end
    end
    
    if params.early_stop == 1 && iter>2
        if iter>i_term
            ss=performance_curve(end-i_term:end)-performance_curve(end);
            if performance_curve(end)~=0
                ss=ss/performance_curve(end);
            end
            if (max(abs(ss))<1e-4)
                const_perf=1;
            end
        end
        
        if exist('accur_max')
            if performance_curve(iter)-accur_max<-1e-4
                first_test=1;
            end
        end
        
        if first_test>(i_term-1) || const_perf
            
            
            if const_perf==0
                iter=iter-(i_term);
            end
            fprintf('Overfitting at iteration %d \n',iter);
            
            Adc=Adc_total{iter};
            X_tr=full(X_total{iter});
            break
        elseif first_test>1
            first_test=first_test+1;
        end
        if performance_curve(iter)>performance_curve(iter-1) && first_test==1
            accur_max=performance_curve(iter-1);
            first_test=first_test+1;
        end
    end
    
    %% ============= Dictionary Update===============
    Adc0=Adc;
    auxD = zeros(size(Adc)) ;
    p=randperm(MP);
    %     p=1:MP;
    thresh_sim=params.sim_thresh;
    for i_c = 1:MP
        c=p(i_c);
        i_nz = abs(X(c,:) ) >  0 ;  % index of non-zero xij
        
        if (isempty(find(i_nz, 1)))
            continue;
        end
        
        E_nz = A(:,i_nz) + Adc(:,c)*X(c,i_nz) ; % compute E_nz by shrinking matrix
        if params.self_rec==0
            id_out=find(i_nz);
        else
            id_out=[];
        end
        
        [ai_up, xj] = ai_update(E_nz, Kyy, X(c,i_nz), Adc, params , H_tr,list_labels_tr,id_out);
        
        % update dictionary column and prune similar atoms
        simi = abs(ai_up'*Kyy*auxD(:,p(1:i_c))) ;
        if find(simi>thresh_sim)
            Adc(:,c)=ai_up;
            X(c,i_nz) = xj;
            
            while(find(simi>thresh_sim))
                [v, i_c0]=max(simi);
                cj=p(i_c0);
                aj=Adc(:,cj);
                i_aj=find(aj);
                if sum(intersect((i_aj),find(i_nz)))==0 || (params.self_rec==1 && sum(intersect((i_aj),find(i_nz))) && numel(i_aj)>1)  % self reconstruction rule
%                 if sum(intersect((i_aj),find(i_nz)))==0 || params.self_rec==1  % self reconstruction rule
                    Adc(:,c)= zeros(n_tr,1) ;
                    i_nzj= i_nz+abs(X(cj,:) ) >  0 ;  % index of non-zero aij
                    
                    X(c,i_nz)=0;
                    A = eye(n_tr) - Adc*X ;
                    E_nzj=A(:,i_nzj) + Adc(:,cj)*X(cj,i_nzj) ; % compute E_nz by shrinking matrix
                    xj=max((aj'*Kyy*E_nzj)/(aj'*Kyy*aj),0);
                    X(cj,i_nzj) = xj;
                    %                     e2=error_kern(Adc,X,Kyy,Kyy,Kyy);
                    %                     e2;
                    break
                else
                    simi(i_c0)=0;
                end
            end
            auxD(:,c)= Adc(:,c);
            
        else
            auxD(:,c) = ai_up ;
            X(c,i_nz) = xj;
        end
        Adc(:,c) = auxD(:,c);
        A = eye(n_tr) - Adc*X ;
        
    end
    %     err_temp=PMSE_kern(Adc,X,Kyy,Kyy,Kyy);
    iter=iter+1;
end
fprintf('Algorithm convereged\n');
% fprintf('=========================\n');
Adc_cvg=Adc;
X_cvg=X;
%% Removing redundancies among the dictioanry elements
disp(' ==================== Pruning the Dictionary =====');
[Adc,err_tr,X]=dic_redunt(Adc_cvg,Kyy,params,X_cvg,0.8);
if params.train_ratio<1
    X_ts=nnkomp_all(Adc,Kyy0,Kzy,Kzz,T0);
    [X_val, err_val]=nnkomp_all(Adc,Kyy0,Khy,Khh,T0);
end
if lcNN
    disp(' ==================== Mapping Accuracy ========');
    [accuracy_tr,conf_mat, i_missed]=LC_classification(X,Adc,H_tr,list_labels_tr);
    fprintf('Accuracy of train data = %3.2f%% \n',accuracy_tr);
    if params.train_ratio<1
        [accuracy_ts,conf_mat_ts, i_missed_ts]=LC_classification(X_ts,Adc,H_tr,list_labels_ts);
        [accuracy_val,conf_mat_val, i_missed_val]=LC_classification(X_val,Adc,H_tr,list_labels_val);
        fprintf('Accuracy of validation data = %3.2f%% \n',accuracy_val);
    end
end
%=== out the results
results.Adc=Adc;
results.kernel=Kyall;
results.Kyy=Kyy0;
if params.train_ratio<1
    results.rec_err=err_val;
else
    results.rec_err=err_tr;
end
results.r_tr=r_tr;
results.X=X;
if params.train_ratio<1
    results.r_ts=r_ts;
    results.r_val=r_val;
    results.X_ts=X_ts;
    results.X_val=X_val;
end
results.list_labels=list_labels;
results.H_tr=H_tr;
if lcNN
    results.accuracy_tr.accuracy=accuracy_tr;
    results.accuracy_tr.conf_mat=conf_mat;
    results.accuracy_tr.i_missed=i_missed;
    if params.train_ratio<1
        results.accuracy_ts.accuracy=accuracy_ts;
        results.accuracy_ts.conf_mat=conf_mat_ts;
        results.accuracy_ts.i_missed=i_missed_ts;
        
        results.accuracy_val.accuracy=accuracy_val;
        results.accuracy_val.conf_mat=conf_mat_val;
        results.accuracy_val.i_missed=i_missed_val;
    end
end