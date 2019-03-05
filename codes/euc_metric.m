function [distance, path1, path2] = euc_metric(x,y,options)
%NONEALIGNMENT Does nothing but a trivial alignment
%   [DISTANCE PATH1 PATH2]=NONEALIGNMENT(X,Y,OPTIONS) produces
%   a trivial alignment in which warping is a simple line.
sx= length(x);
sy= length(y);
if sx ~= sy
    sm= min(length(x),length(y));
    ix=ceil(linspace(1,sx,sm));ix(end)=sx;
    iy=ceil(linspace(1,sy,sm));iy(end)=sy;
    x=x(ix);
    y=y(iy);
end
distance = sqrt(sum((x - y).^2));
path1 = 1:length(x);
path2 = 1:length(y);
end