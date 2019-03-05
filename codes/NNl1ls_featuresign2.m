function [Xout,fobj] = NNl1ls_featuresign2 (H, Cv, gamma, qp_solve,T0,idy, alfasum, Xinit)
% The feature-sign search algorithm
% L1-regularized least squares problem solver
%
% This code solves the following problem:
%
%    minimize_s 0.5*||y - A*x||^2 + gamma*||x||_1
%    minimize_s 0.5*x'*H*x + c'*x + gamma*||x||_1
%
% The detail of the algorithm is described in the following paper:
% 'Efficient Sparse Codig Algorithms', Honglak Lee, Alexis Battle, Rajat Raina, Andrew Y. Ng,
% Advances in Neural Information Processing Systems (NIPS) 19, 2007
%
% Written by Honglak Lee <hllee@cs.stanford.edu>
% Copyright 2007 by Honglak Lee, Alexis Battle, Rajat Raina, and Andrew Y. Ng

warning('off', 'MATLAB:divideByZero');

use_Xinit= false;
if exist('Xinit', 'var')
    use_Xinit= true;
end
if ~exist('alfasum')
    alfasum=0;
end
if ~exist('qp_solve')
    qp_solve=0;
end;

if sum(Cv==0) && numel(Cv)==1
    Cv=zeros(size(H,1),1);
end

if ~exist('T0')
    T0=length(Cv);
end;

if ~exist('idy')
    idy=[];
end
% if T0<length(Cv)
%     gamma=0;
% end
% Xout= zeros(size(A,2), size(Y,2));
AtA = H; % =H
AtY = Cv; % =c'

% rankA = rank(AtA);


% for i=1:size(Y,2)
%     if mod(i, 100)==0, fprintf('.'); end %fprintf(1, 'l1ls_featuresign: %d/%d\r', i, size(Y,2)); end

if use_Xinit
    %         rankA = min(size(Xinit,1)-10, size(A,2)-10);
    T0=length(Cv)+1;
    rankA=size(Xinit,1);
    idx1 = find(Xinit~=0);
    maxn = min(length(idx1), rankA);
    xinit = zeros(size(Xinit));
    xinit(idx1(1:maxn)) =  Xinit(idx1(1:maxn));
    [Xout, fobj]= NNls_featuresign_sub (AtA, AtY, gamma, qp_solve, T0,idy, alfasum, xinit);
else
    [Xout, fobj]= NNls_featuresign_sub (AtA, AtY, gamma, qp_solve, T0,idy,alfasum);
end
% end
% fprintf(1, '\n');

warning('on', 'MATLAB:divideByZero');

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [x, fobj] = NNls_featuresign_sub (H, Cv, gamma, qp_solve, T0,idy, alfasum, xinit)

[L,M] = size(H);

rankA = M;% ?? not sure!!

% Step 1: Initialize
usexinit = false;
if ~exist('xinit', 'var') || isempty(xinit)
    xinit= [];
    x= sparse(zeros(M,1));
    theta= sparse(zeros(M,1));
    act= sparse(zeros(M,1));
    allowZero = false;
else
    % xinit = [];
    x= sparse(xinit);
    theta= sparse(sign(x));
    act= sparse(abs(theta));
    usexinit = true;
    allowZero = true;
end

fname_debug = sprintf('../tmp/fsdebug_%x.mat', datestr(now, 30));

fobj = 0; %fobj_featuresign(x, A, y, AtA, Aty, gamma);

ITERMAX=1000;
optimality1=false;
optimality2=0;
for iter=1:ITERMAX
    % check optimality0
    if optimality2==1 || numel(find(x)) >= T0
        break;
    end
    if iter>20
        iter;
    end
    act0=act;
    act0(idy)=1;
    act_indx0 = find(act0 == 0);
    %     if ~isempty(idy)
    %         if find(act_indx0==idy)
    %             act_indx0(find(act_indx0==idy))=[];
    %         end
    %     end
    grad = H*sparse(x) + Cv;% change H*x+c
    theta = sign(x);
    
    %     optimality0= false;
    % Step 2
    %     while ~isempty(act_indx0)
    %     [mx,indx] = max (abs(grad(act_indx0)));
    [mn,indx] = min (grad(act_indx0));
    %     [mn,indx] = max (grad(act_indx0));
    if ~isempty(mn)
        if max(grad(act_indx0))>0 && sum(x)==0
            mn;
        end
    end
    
    if ~isempty(mn) && (mn <= -gamma) && (iter>1 || ~usexinit)
        %     if ~isempty(mn) && (mn >= 0) && (iter>1 || ~usexinit)
        act(act_indx0(indx)) = 1;
        theta(act_indx0(indx)) = 1;
        usexinit= false;
    else
        %         optimality0= true; %dummy
        save temp_gi0
        error('gi==0 problem! it is saved in temp_gi0')
    end
    act_indx1 = find(act == 1);
    
    if length(act_indx1)>rankA
        warning('sparsity penalty is too small: too many coefficients are activated');
        fobj = fobj_featuresign(x,H, Cv, gamma);
        return;
    end
    
    if isempty(act_indx1) %length(act_indx1)==0
        % if ~assert(max(abs(x))==0), save(fname_debug, 'A', 'y', 'gamma', 'xinit'); error('error'); end
        if allowZero, allowZero= false; continue, end
        return;
    end
    
    % Step 3: feature-sign step
    % if ~assert(length(act_indx1) == length(find(act==1))), save(fname_debug, 'A', 'y', 'gamma', 'xinit'); error('error'); end
    k=0;
    %     while 1
    %         k=k+1;
    if k> 100
        %         k;
    end
    
    if k>ITERMAX
        warning('Maximum number of iteration reached. The solution may not be optimal');
        fobj = fobj_featuresign(x,H, Cv, gamma);
        % save(fname_debug, 'A', 'y', 'gamma', 'xinit');
        %             return;
        break
    end
    
    if isempty(act_indx1) % length(act_indx1)==0
        % if ~assert(max(abs(x))==0), save(fname_debug, 'A', 'y', 'gamma', 'xinit'); error('error'); end
        if allowZero, allowZero= false; break, end
        return;
    end
    
    % Step 3: feature-sign step
    [x, theta, act, act_indx1, optimality1, lsearch, fobj, optimality2] = compute_FS_step (x, H, Cv, theta, act, act_indx1, gamma, qp_solve, alfasum);
    
    % Step 4: check optimality condition 1
    %         if optimality1
    grad = H*sparse(x) + Cv;% change H*x+c
    if abs(grad(find(x))) >1e-4
        %                 optimality1;
        %             else
        error('g(x~=0) is not 0')
    end
    %             if abs(grad(find(x==0)))<=(1.01)*gamma
    if sum(grad(find(x==0))+1e-4 <0)==0
        optimality2=1;
    end
    %             break;
    %             end
    %         end;
    if lsearch >0  % dummy
        %             continue;
    end;
    
    %     end
end

if iter >= ITERMAX
    warning('Maximum number of iteration reached. The solution may not be optimal');
    % save(fname_debug, 'A', 'y', 'gamma', 'xinit');
end

if 0  % check if optimality
    act_indx1 = find(act==1);
    grad = AtA*sparse(x) - Aty;
    norm(grad(act_indx1) + gamma.*sign(x(act_indx1)),'inf')
    find(abs(grad(setdiff(1:M, act_indx1)))>gamma)
end

fobj = fobj_featuresign(x,H, Cv, gamma);

return;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, theta, act, act_indx1, optimality1, lsearch, fobj, optimality2] = compute_FS_step (x, H, Cv, theta, act, act_indx1, gamma, qp_solve, alfasum)
x2 = x(act_indx1);
% A2 = A(:, act_indx1);
theta2 = theta(act_indx1);
H2 = H(act_indx1, act_indx1);
Cv2=(Cv(act_indx1) + gamma.*theta2);
% Cv2=(Cv(act_indx1) + gamma.*theta2/2);
obj=@(x) 0.5*x'*H2*x+Cv2'*x;
obj0=obj(x2);

% if f_NN==0
% call matlab optimization solver..
% x_new = - H2 \ Cv2; % change to NNQP
x_new = - pinv(H2) * Cv2; % change to NNQP
eh=eig(H2);
if abs(min(eh)/max(eh))*100>.1 && min(eh)<0
    disp('non-neg eigens for h2');
end
if sum(x_new)==0 && qp_solve==0
    %     save temp_x_zero
    display('sum(x_new)==0')
end
% eg=eig(H2);
% if find(abs(eg)<1e-4)
%     H2;
% end
% x_new2 = sparse(NNQP_SC(H2,Cv2,1e-5));
%     x_new = H2 \ ( - Cv(act_indx1) - gamma.*theta2 ); % change to NNQP
% f_obj1=0.5*x_new'*H2*x_new+( Cv(act_indx1) + gamma.*theta2/2 )'*x_new;%+0.5*y'*y;
% else
%     Cv2=(Cv(act_indx1) + gamma.*theta2/2);
options_fmin = optimoptions('fmincon','Algorithm','sqp','Display','off');

if qp_solve==1 || alfasum
    %         x_new = sparse(NNQP_SC(H2,Cv2,1e-5));
    %     else
    %     options_qp = optimoptions('quadprog','Algorithm','active-set','Display','off');
    %         options = optimoptions('quadprog','Algorithm','trust-region-reflective','Display','off');
    %         x_new = quadprog(H2,Cv2,[],[],[],[],0*ones(size(H2,1),1),[],[],options_qp);
    if alfasum
        x_new = fmincon(@(x) 0.5*x'*H2*x+Cv2'*x,rand(size(Cv2)),[],[],ones(size(Cv2))',1,0*ones(size(H2,1),1),[],[],options_fmin);
    else
        x_new = fmincon(@(x) 0.5*x'*H2*x+Cv2'*x,rand(size(Cv2)),[],[],[],[],0*ones(size(H2,1),1),[],[],options_fmin);
    end
    
end

% [x_new2,f_obj2,i_end] = NNQP_SC(H2,Cv2,1e-5);
% f_obj2=0.5*x_new2'*H2*x_new2+(Cv(act_indx1) + gamma.*theta2/2 )'*x_new2;%+0.5*y'*y;
% opts.POSDEF=true; opts.SYM=true; % RR
% x_new = linsolve(AtA2, ( Aty(act_indx1) - gamma.*theta2 ), opts); % RR
% if find(x_new<0)
%     x_new;
% end
% if f_obj2>f_obj1
%     f_obj2;
% end
% x2f=full(x2);
optimality1= false;  %dummy
optimality2= 0;
if (sign(x_new) >= 0) %&& sum(x2~=0)
    if size(x_new,1)>1
        x_new;
    end
    optimality1= true;
    x(act_indx1) = x_new;
    fobj = obj(x_new); %fobj_featuresign(x, A, y, AtA, Aty, gamma);
    lsearch = 1;
    return;
else  % projection
    if sum(x_new<0)>1
        progress = (0 - x2)./(x_new - x2);
        progress2 = (0 - x_new)./(x2-x_new);
%         progress(find(x_new<0))
%         progress2(find(x_new<0))
%         disp('sum(x_new<0)>1')
    end
    xx=fmincon(@(x) 0.5*x'*H2*x+Cv2'*x,rand(size(Cv2)),[],[],[],[],0*ones(size(H2,1),1),[],[],options_fmin);
    objxx=obj(xx);
    i1=find(xx==0);
    i2=find(x_new<0);
    i_p=0;
    if numel(i2)~=numel(i1)
        i_p=1;
    elseif sum(i1==i2) ~= numel(i2)
        i_p=1;
    end
    
%     if  i_p==1
%         grad=H2*x_new + Cv2;
%         [val , ig]=sort(grad);
%         ig;
% %         disp('sum(i1==i2) ~= numel(i2)')
%         x_new=xx;
%         i_p=0;
%     end
    while (1)
        if isempty(act_indx1)
            error('isempty(act_indx1)')
        end
        progress = (0 - x2)./(x_new - x2);
        i_n=find(x_new<0);
        p_neg=progress(find(x_new<0));
        [v i_p]=min(p_neg);
        
        remove_idx=i_n(i_p);
        %         remove_idx=find(x_new<=0);
        
        x(act_indx1(remove_idx))=0;
        theta(act_indx1(remove_idx))=0;
        act(act_indx1(remove_idx))=0;
        act_indx1(remove_idx)=[];
        
        
        theta2 = theta(act_indx1);
        H2 = H(act_indx1, act_indx1);
        Cv2=(Cv(act_indx1) + gamma.*theta2);
        % Cv2=(Cv(act_indx1) + gamma.*theta2/2);
        obj=@(x) 0.5*x'*H2*x+Cv2'*x;
        
        x3 = - pinv(H2) * Cv2; % change to NNQP
        if sum(x3<0)==0
            x_new=x3;
            obj2=obj(x3);
            if  obj2>obj0 
                disp('not optimal x_new')
            end
            if numel(find(xx>0))~=numel(x_new)
                disp('different x_new')
            elseif norm(xx(xx>0)-x_new) > 0.1 *norm(x_new) || obj2>obj0
                disp('not optimal x_new')
            end
            
            
            x(act_indx1) = x_new;
            fobj = obj(x_new);
            lsearch=1;
            optimality1=1;
            return;
        else
            x2 = x(act_indx1);
            x_new=x3;
%             disp('negative x3')
        end
    end
    
end

% do line search: x -> x_new
progress = (0 - x2)./(x_new - x2);
lsearch=0;
%a= 0.5*sum((A2*(x_new- x2)).^2);
a= 0.5*(x_new- x2)'*H2*(x_new- x2);
b= (x2'*H2*(x_new- x2) + (x_new- x2)'*Cv(act_indx1));
fobj_lsearch = gamma*sum(abs(x2));
[sort_lsearch, ix_lsearch] = sort([progress',1]);
remove_idx=[];
for i = 1:length(sort_lsearch)
    t = sort_lsearch(i);
    if t<=0 | t>1
        continue;
    end
    s_temp= x2+ (x_new- x2).*t;
    fobj_temp = a*t^2 + b*t + gamma*sum(abs(s_temp));
    if fobj_temp < fobj_lsearch
        fobj_lsearch = fobj_temp;
        lsearch = t;
        if t<1  remove_idx = [remove_idx ix_lsearch(i)]; end % remove_idx can be more than two..
    elseif fobj_temp > fobj_lsearch
        break;
    else
        if (sum(x2==0)) == 0
            lsearch = t;
            fobj_lsearch = fobj_temp;
            if t<1  remove_idx = [remove_idx ix_lsearch(i)]; end % remove_idx can be more than two..
        end
    end
end

% if ~assert(lsearch >=0 && lsearch <=1), save(fname_debug, 'A', 'y', 'gamma', 'xinit'); error('error'); end

if lsearch >0
    % update x
    x_new = x2 + (x_new - x2).*lsearch;
    x(act_indx1) = x_new;
    theta(act_indx1) = sign(x_new);  % this is not clear...
end

% if x encounters zero along the line search, then remove it from
% active set
% if lsearch<1 & lsearch >0 %???
if lsearch<1
    %remove_idx = find(x(act_indx1)==0);
    remove_idx = find(abs(x(act_indx1)) < eps);
    x(act_indx1(remove_idx))=0;
    
    theta(act_indx1(remove_idx))=0;
    act(act_indx1(remove_idx))=0;
    act_indx1(remove_idx)=[];
end

%--- optimality
% grad = H*sparse(x) + Cv;% change H*x+c
% if norm(grad(find(x))-gamma*sign(x(find(x))))/norm(grad(find(x)))*100 <1
%     optimality1= true;
%     x(act_indx1) = x_new;
%     fobj = 0; %fobj_featuresign(x, A, y, AtA, Aty, gamma);
%     lsearch = 1;
%     return
% end

% if (sign(x_new) == full(sign(x2)))
%     if size(x_new,1)>1
%         x_new;
%     end
%     optimality1= true;
%     x(act_indx1) = x_new;
%     fobj = 0; %fobj_featuresign(x, A, y, AtA, Aty, gamma);
%     lsearch = 1;
%     return;
% end


fobj_new = 0; %fobj_featuresign(x, A, y, AtA, Aty, gamma);

fobj = fobj_new;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f, g] = fobj_featuresign(x, H, Cv, gamma)

% f= 0.5*norm(y-A*x)^2;
f= 0.5*x'*H*x+Cv'*x;
f= f+ gamma*norm(x,1);

if nargout >1
    g= H*x + Cv;
    g= g+ gamma*sign(x);
end

return;

%%%%%%%%%%%%%%%%%%%%%

function retval = assert(expr)
retval = true;
if ~expr
    % error('Assertion failed');
    warning ('Assertion failed');
    retval = false;
end
return
