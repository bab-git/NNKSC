function [Adc, Adc_nrm]=dic_clean(Adc,thresh)
%Pruning elements of Adc based on a given threshold
if ~exist('thresh')
    thresh=0.25;
end
Adc_nrm=Adc./repmat(max(abs(Adc)),[size(Adc,1) 1]);
Adc_nrm(isnan(Adc_nrm))=0;
Adc(abs(Adc_nrm)<thresh)=0;