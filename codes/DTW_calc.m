function [alldistances,alldistances_feats]=DTW_calc(raw_data,Alignment,DTW_nrm,w_win)
% Computing the pair-wise alignment between motion data to obtain distance
% matrix describing the dataset

% INPUTS:
  % raw_data: raw data in the form of raw_data{i}.data
  % Alignment: DTW=DTW alighnment , Euc=Euclidean distance calculation
  % DTW_nrm:  1=normalizing DTW output according to the length of the path
  % w_win:    warping window for DTW  chosen between (0 1] and multiplied by the signal's length
  
% OUTPUTS: 
    % alldistances= Computed pair-wise distances between data points.
    % alldistances_feats= Computed the pair-wise distances for each dimension
    % of the data individually

if ~exist('Alignment')
    Alignment='DTW';
end

if ~exist('w_win')
    w_win=1;
end

if ~exist('DTW_nrm')
    DTW_nrm=1;
end

normalize=1;  % normalizing input signals
sub_samp=1;   %subsampleing of time series if required
%% checking the data format
if ~isfield(raw_data{1},'data')
    error('Please put each motion sequence "i" in form of "raw_data{i}.data"')
end

n_data=length(raw_data);
if size(raw_data{1}.data,1)==size(raw_data{2}.data,1) && size(raw_data{3}.data,1)==size(raw_data{2}.data,1)
    feat_dim=1;
elseif size(raw_data{1}.data,2)==size(raw_data{2}.data,2) && size(raw_data{3}.data,2)==size(raw_data{2}.data,2)
    feat_dim=2;
else
    error('each "dataset.data{i}" should contain a multi-dimension time-series')
end
nf=size(raw_data{1}.data,feat_dim);

alldistances=zeros(n_data,n_data,nf);
%%
for is=1:n_data
    
    
    for js=is+1:n_data
        for ks=1:nf
            if feat_dim==1
                x=raw_data{is}.data(ks,:);
                y=raw_data{js}.data(ks,:);
            else
                x=raw_data{is}.data(:,ks);
                y=raw_data{js}.data(:,ks);
            end
            if sub_samp~=1
                x=resample(x,1,1/sub_samp);
                y=resample(y,1,1/sub_samp);
            end
            if normalize==1
                x=zscore(x);
                y=zscore(y);
            end
            %% =================== Alignment                       
            if strcmpi(Alignment,'Euc')
                [alldistances(is,js,ks), path1, path2] = euc_metric(x,y,[]);
            elseif strcmpi(Alignment,'DTW')
                
                win=ceil(max(length(x),length(y))*w_win); % window
                
                
                [dist,~,path_l]=dtw_freq(x',y',win);
                
                if DTW_nrm
                    dist=dist/path_l;
                end
                alldistances(is,js,ks)=dist;
                
                
            end
        end
    end
    fprintf('%3.2f%% of data is processed \n ',is/n_data*100);  %if (mod(iter,30)==0), fprintf('\n'); end
end
alldistances_feats=dist_correct2(alldistances);
alldistances=sum(alldistances_feats,3);
