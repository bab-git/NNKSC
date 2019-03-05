function [A, x, MSE] = ai_update( E, Ky, initx, Ad,params, C,labels,id_out,loop_drw)
% This function updates ai (a column of the dictionay matrix) and xj (the entries in X which use ai)

[d, n] = size( E );

numIter=params.loop_a;

if ~exist('initx','var')
    initx = rand( 1, n );  % init A randomly
end

exitTol = 0.1;

if ~exist('loop_drw','var')
    loop_drw = 0;
end

x = initx;
XtX=Ky;
if ~exist('id_out')
    id_out=[];
end
xd=diag(XtX);
if params.class_a   % forcing a to use same class of data
    numIter=100;
end
iter=1;
%%
while iter <= numIter
    if loop_drw>0
        if rem(iter,10)==0
            fprintf( '\n  Mean Squared Error: %d \n', MSE(iter) );
            fprintf( '\n  step %d', iter );
        else
            fprintf( ' %d', iter );
        end
    end
    
    [A, fobj]= NAMP_A(Ad, x, Ky, XtX, xd, E,params,id_out);
    
    if params.class_a
        c_v=C*A;
        [val ind]=max(c_v);
        if numel(find(c_v))==1
            numIter=iter;
        else
            id_out=[id_out find(labels~=ind)];
        end
    end
    
    x=max((A'*Ky*E)/(A'*Ky*A),0);
    norm_k=error_E(E,Ky,A,x);
    
    MSE(iter)=norm_k;
    if MSE(iter)< exitTol
        break
    end
    if iter>5
        if max(abs(MSE(end-5:end)/MSE(end)-1))<exitTol
            break
        end
        
    end
    iter=iter+1;
end
%%
if sum(A)>0
    norm_phia=sqrt(A'*Ky*A);
    x=x*norm_phia;
    A=A/norm_phia;
end
A=sparse(A);
if loop_drw==1
    fprintf( '\n\n  finished after %d iterations.' ,iter );
    fprintf( '\n  Mean Squared Error: %d \n', MSE(end) );
end
