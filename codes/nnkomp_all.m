function [X,res_all]=nnkomp_all(Adc,Kyy,Kzy,Kzz,T0,SR)
tr_data=0;
if ~exist('SR')
    SR=1;
end
if size(Kzy,1)==size(Kzy,2)
    if norm(Kzy-Kzy')<1e-8
        tr_data=1;  % it's reconstruction of the training data
    end
end

Adc_temp=Adc;
G0=Adc'*Kyy*Adc;
B0=Adc'*Kzy';
G=G0;
B=B0;
for i_y=1:size(Kzz,1)
    
    if tr_data
        %-- removing direct solutions xi=yi
        Adc_temp=Adc;
        id=find(Adc(i_y,:));        
        id2=sum(Adc(:,id)>0);
        id3=id(find(id2==1));
        if SR==0   % no self reconstruction
            id3=[id3 id];
        end
        Adc_temp(:,id3)=0;
        G=G0;G(id3,:)=0;G(:,id3)=0;
        B=B0;B(id3,:)=0;
    end
    
    [x, res_x]= nnKomp(Adc_temp, Kyy, Kzy(i_y,:),Kzz(i_y,i_y),T0,G,B(:,i_y));
    X(:,i_y)=x;
    res_all(i_y)=res_x;
end