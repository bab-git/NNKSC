function [fig, plotHandles] = plotData_MP2(Y,r_tr, labels_class, plotAxis, prot_y,show_num,offset)
% pass figHand = 0 or don't pass it if a new figure should be opened.
X=Y(r_tr,:);
labels_tr=labels_class(r_tr);
% old call:
% function [fig, plotHandles] = plotData(X, labels, idx, figHand, axisOff)
off_rate=1;
if ~exist('plotAxis')
    plotAxis = 0;
elseif isempty(plotAxis)
    plotAxis = 0;
end

if ~exist('show_num')
    show_num=0;
end


if ~exist('offset')
    offset = 0;
elseif isempty(offset)
    offset = 0;
end

% if (nargin >= 7) && (~isempty(idx))
%     labels_MP = labels_MP(idx);
%     X = X(idx,:);
% end

% figures for paper?
forPaper = 0;

if forPaper
    markerSize   = 10;
    fontSizeFigs = 28;
else
    markerSize   = 20;
    fontSizeFigs = 10;
end


% transpose the data if necessary
if size(X,2) > 3
    if size(X,1) <= 3
        warning('Transposing the data matrix to have the instances in the rows');
        X = X';
    else
        error('Data should have maximaly 3 features.');
    end
end

% compute number of classes
% uniqueMP = unique(labels_MP);
nMP     = length(prot_y);
nData   = size(Y,1);
nCl     = max(labels_tr);
dim          = size(X,2);


normalLabels = 1;
if (length(unique(labels_tr)) > 300) || (sum(abs(round(labels_tr) - labels_tr)) > 0)
    % Assuming a regression and plotting continuois!
    origLabels   = labels_tr;
    normalLabels = 0;
    warning('Assuming a regression and plotting continuous!');
    %scaleTo = max(1000, length(labels));
    scaleTo = 1000;
    labels_MP  = labels_MP - min(labels_MP);
    scale   = scaleTo/max(labels_MP);
    labels_MP  = round(labels_MP*scale) + 1;
    
    uniqueMP = unique(labels_MP);
    nMP     = scaleTo+1;
end




if normalLabels
    % '.', 'x', '+', are left out
    symbs = {'o', 's', 'd', 'v', '^', '<', '>', 'p', 'h', '*'};
    %symbs = {'o', '*', 'd', 'v', '^', '<', '>', 'p', 'h'}; % tmp
    if nCl > length(symbs)
        % if more classes then symbols exists, repeat the latter
        symbs = repmat(symbs, 1, ceil(nCl/length(symbs)));
    end
    %symbs = {'.', 'o'};
else
    % use only one marker
    symbs = repmat({'o'}, nCl, 1);
end

if normalLabels
    cm = hsv(nCl);
    if (nCl == 2)
        cm = [1 0 0; 0 0 1];
    end
else
    cm = jet(nCl);
end


% if (nargin >= 8) && (~isempty(figHand))
if exist('figHand')
    %         fig = myFigure;
    %     else
    fig = figHand;
    %     end
else
    fig = figure;
end

hold on; box on;
% if axisOff
%     axis off;
% end

set(fig, 'color', 'white');

legendStr = cell(nCl,1);
for i=1:nCl
    legendStr{i} = ['Cls ' int2str(i)];
end




% save handles
% plotHandles = zeros(nMP,1);
% if normalLabels
% used_data=[];
% used_data=zeros(nData,1);
% i_hand=0;
i_cls_did=[];
for iy=1:length(prot_y)
    i_mp=prot_y{iy};
    if isempty(i_mp)
        continue
    end
    %     for i_y=1:numel(temp_i)
    %         iy=temp_i(i_y);
    i_cls=labels_tr(iy);
    face_clr=cm(i_cls,:);
    tmp = plot(X(iy,1), X(iy,2), 'o','MarkerSize', markerSize,'MarkerEdgeColor',cm(i_cls,:));
    text(  X(iy,1),  X(iy,2)+offset, num2str(i_mp), 'Color', cm(i_cls,:), 'FontWeight', 'normal', 'FontSize', fontSizeFigs, 'HorizontalAlignment', 'center','fontweight','bold');
    
    if show_num
        text(  X(iy,1),  X(iy,2)-off_rate, num2str(iy), 'Color', [0 0 0], 'FontWeight', 'normal', 'FontSize', fontSizeFigs, 'HorizontalAlignment', 'center');
    end
    
    if sum(i_cls==i_cls_did)==0
        i_cls_did=[i_cls_did i_cls];
        %             i_hand=i_hand+1;
        plotHandles(i_cls) = tmp;
    end
    
end

% i_hand=0;
% i_cls_did=0;
for iy=1:nData
    i_cls=labels_class(iy);
    if i_cls==0
        continue
    end
    i_used=0;
    if sum(r_tr==iy)
        i_us=find(r_tr==iy);
        if ~isempty(prot_y{i_us})
            i_used=1;
        end
        
    end
    if i_used==0
        
        tmp = plot(Y(iy,1), Y(iy,2), symbs{i_cls},'MarkerSize', 10,'MarkerFaceColor',cm(i_cls,:),'MarkerEdgeColor',cm(i_cls,:));
        
        %         if (i_cls_did+1==i_cls)
        %             i_cls_did=i_cls_did+1;
        %             i_hand=i_hand+1;
        %             plotHandles(i_hand) = tmp;
        %         end
    end
end

if plotAxis
    grid on;
    xlabel('Dim 1',  'FontSize', fontSizeFigs, 'FontWeight', 'bold'); % 22
    ylabel('Dim 2',  'FontSize', fontSizeFigs, 'FontWeight', 'bold');
    if dim == 3
        zlabel('Dim 3',  'FontSize', fontSizeFigs, 'FontWeight', 'bold');
    end
end


% determine the string for the labeling automatically
legendStr = cell(nCl,1);
for i=1:nCl
    legendStr{i} = ['Cls ' int2str(i)];
end

% if normalLabels
i_c=unique(i_cls_did);
legendStr=legendStr(i_c);
plotHandles=plotHandles(i_c);
h_legend = legend(plotHandles, legendStr);
set(h_legend,'FontSize',fontSizeFigs); set(h_legend,'location','best');
% end

if ~plotAxis
    if dim==2
        set(gca,'XTickLabel','','yTickLabel','','xtick',-100:10:-110,'ytick',-100:10:-110);
    elseif dim==3
        set(gca,'XTickLabel','','yTickLabel','', 'zTickLabel', '','xtick',-100:10:-110,'ytick',-100:10:-110,'ztick',-100:10:-110);
    end
end

set(gca, 'FontSize', fontSizeFigs, 'FontWeight', 'bold' );

if ~normalLabels
    colorbar;
    cbar_handle = findobj(fig,'tag','Colorbar');
    set(cbar_handle, 'YTick', []);
end

if show_num
    descr = {'Plot shows how the dictionary primitives are constructed from the data:';
        'Colors --> classes of data';
        '                            Colored Numbers --> Dictionary Primitives';
        '            Black Numbers --> Data samples';
        };
else
    descr = {'Plot shows how the dictionary prmitives are constructed from the data:';
        'Colors --> classes of data';
        '                            Colored Numbers --> Dictionary Primitives';    
        };
%       '                                   Colored dots --> Unused data in the dictionary';
end
title(descr);
% set(gcf,'title','Myfig')
% set view of an old figure
% [AL, EL] = view; view(AL, ELL);
% currAx = axis; axis(currAx);


