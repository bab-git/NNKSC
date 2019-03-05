% maps the prototypes and data together
function [fig, plotHandles] = plotData_MP_mar(X, labels,MP, plotAxis,show_num,font,markerSize, labels_all, fig_Handle)
% pass fig_Handle = 0 or don't pass it if a new figure should be opened.

off_rate=0;

if ~exist('fig_Handle')
    fig_Handle=0;
end

if ~exist('show_num')
    show_num=0;
end

if ~exist('plotAxis') || isempty(plotAxis)
    plotAxis = 0;
end
n_data=length(labels)-MP;
% if (nargin >= 6) && (~isempty(x_ind))
%     labels = labels(x_ind);
%     X = X(x_ind,:);
% end

% figures for paper?
% forPaper = 0;

% if forPaper
%     markerSize   = 10;
%     fontSizeFigs = 28;
% else
if ~exist('markerSize') || isempty(markerSize)
    markerSize   = 10;
end
%     fontSizeFigs = 7;
% end
if exist('font') && ~isempty(font)
    %         markerSize   = font;
    fontSizeFigs = font-5;
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
uniqueLabels = unique(labels);
uniqueLabels2 = unique(labels_all);
nClasses     = numel(uniqueLabels);
nClasses2     = numel(uniqueLabels2);
dim          = size(X,2);

normalLabels = 1;
% if 0*(length(unique(labels)) > 30) || (sum(abs(round(labels) - labels)) > 0)
%     % Assuming a regression and plotting continuois!
%     origLabels   = labels;
%     normalLabels = 0;
%     warning('Assuming a regression and plotting continuous!');
%     %scaleTo = max(1000, length(labels));
%     scaleTo = 1000;
%     labels  = labels - min(labels);
%     scale   = scaleTo/max(labels);
%     labels  = round(labels*scale) + 1;
%
%     uniqueLabels = unique(labels);
%     nClasses     = scaleTo+1;
% end




if normalLabels
    % '.', 'x', '+', are left out
%     symbs = {'o', 's', 'd', 'v', '^', '<', '>', 'p', 'h', 'x','+'};
%     symbs = {'o', '*', 'd', 'v', '^', '<', '>', 'p', 'h'}; % tmp
    symbs = {'s','^','o', '>', 'd', 'v', 'p', '<', 'h'}; % tmp
    if 0*(nClasses2 > length(symbs))
        % if more classes then symbols exists, repeat the latter
        %         symbs = repmat(symbs, 1, ceil(nClasses/length(symbs)));
        symbs = repmat(symbs, 1, ceil(nClasses2/length(symbs)));
    end
    %symbs = {'.', 'o'};
else
    % use only one marker
    symbs = repmat({'o'}, nClasses, 1);
end

% if normalLabels
if nClasses>6
    cm = hsv(nClasses);
    if (nClasses == 2)
        cm = [1 0 0; 0 0 1];
    end
else
%     cm={'r','g','k','b','y','w'};
    cm=[1 0 0;1 1 0;0 1 0;0 1 1;0 0 1;1 1 1];
end
% else
%     cm = jet(nClasses);
% end


if fig_Handle~=0
    fig = fig_Handle;
else
    fig = figure;
end

hold on; box on;
% if axisOff
%     axis off;
% end

set(fig, 'color', 'white');

labels_data=labels(1:n_data);
% save handles
plotHandles = zeros(nClasses2,1);
i_leg=0;
leg_ind=[];
if normalLabels
    for i=1:nClasses
        %         if dim == 2
%         symbs=symbs(randperm(length(symbs)));
        iy_all=find(labels_data==i);
        labels_2=unique(labels_all(iy_all));
        for i2=1:length(labels_2)
            i_leg=i_leg+1;
            i_obj=labels_2(i2);
            color_ind=i;
            iy=find(labels_all==i_obj);
%             if length(unique(labels_all(iy)))>1
%                 error('if length(unique(labels_all(iy)))>1')
%             end
            tmp = plot(X(iy,1), ...
                X(iy,2), ...
                symbs{i2}, 'MarkerSize', markerSize, 'MarkerFaceColor', cm(color_ind,:), ...
                'Color', cm(color_ind,:),'MarkerEdgeColor', 'k','LineWidth',2);
            if show_num
                for i_num=1:numel(iy)
                    text(  X(iy(i_num),1),  X(iy(i_num),2)-off_rate, num2str(iy(i_num)), 'Color', [0 0 0], 'FontWeight', 'bold', 'FontSize', fontSizeFigs, 'HorizontalAlignment', 'center');
                end
            end
            if ~isempty(tmp)
                leg_ind=[leg_ind i_obj];
                plotHandles(i_leg) = tmp;
            end
        end
        iy=n_data+find(labels(end-MP+1:end)==color_ind);
        plot(X(iy,1), ...
            X(iy,2), ...
            'p', 'MarkerSize', markerSize+8, 'MarkerEdgeColor', 'k','MarkerFaceColor', cm(color_ind,:),'LineWidth',3);
        if show_num
            i_prot=find(labels(end-MP+1:end)==color_ind);
            for i_num=1:numel(iy)
                text(  X(iy(i_num),1),  X(iy(i_num),2)-off_rate, num2str(i_prot(i_num)), 'Color', [0 0 0], 'FontWeight', 'bold', 'FontSize', fontSizeFigs-4, 'HorizontalAlignment', 'center');
            end
        end
        
        %             plot(X(iy,1), ...
        %                 X(iy,2), ...
        %                 symbs{i}, 'MarkerSize', markerSize+2, 'MarkerEdgeColor', 'k');
        %             plot(X(iy,1), ...
        %                 X(iy,2), ...
        %                 symbs{i}, 'MarkerSize', markerSize+4, 'MarkerEdgeColor', 'k');
        
        
        %         elseif dim == 3
        %             tmp = plot3(X(labels==uniqueLabels(i),1), ...
        %                 X(labels==uniqueLabels(i),2), X(labels==uniqueLabels(i),3), ...
        %                 symbs{i}, 'MarkerSize', markerSize, 'MarkerFaceColor', cm(i,:), ... % paper version
        %                 'Color', cm(i,:));
        %             if ~isempty(tmp)
        %                 plotHandles(i) = tmp;
        %             end
        %         end
    end
else
    for i=1:numel(uniqueLabels)
        tmp = plot(X(labels==uniqueLabels(i),1), ...
            X(labels==uniqueLabels(i),2), ...
            symbs{i}, 'MarkerSize', markerSize, 'MarkerFaceColor', cm(uniqueLabels(i),:), ...
            'Color', cm(uniqueLabels(i),:));
        if ~isempty(tmp)
            plotHandles(i) = tmp;
        end
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


%range = [-60 60 -70 50];
%axis(range);
% or as an alternative argument for lengend:
%  legend(h, strcat('contig_',strtrim(cellstr(num2str(unique(Labels'))))), 'Location', 'BestOutside', 'Interpreter', 'none')
%  where h is a series of handles returned by plot: h(i) = plot(...).
%h_legend = legend('1', '2', '3', '4', '5', '6', '7', '8', '9', '0'); set(h_legend,'FontSize',16);
%h_legend = legend('ACA', 'ACC');
%h_legend = legend('intron', 'exon');
%h_legend = legend('C1', 'C2', 'C3');
%h_legend = legend(plotHandles, {'C1', 'C2', 'C3'});
% pos = [0.8818 0.0319 0.0840 0.5177];
%pos = [0.0316 0.0788 0.0970 0.6708];
%set(h_legend, 'Position', pos);
%range = axis;
%axis(range);

% legend handle bekommen
% tmp = findobj(gcf,'Type','axes','Tag','legend')


% determine the string for the labeling automatically
legendStr = cell(nClasses2,1);
for i=1:nClasses2
%     legendStr{i} = ['Cls ' int2str(i)];
    legendStr{i} = ['Cls ' int2str(leg_ind(i))];
end
%legendStr = {'ACA', 'ACC'};

if normalLabels
    h_legend = legend(plotHandles, legendStr);
    set(h_legend,'FontSize',fontSizeFigs+5); set(h_legend,'location','best');
end

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

% set view of an old figure
% [AL, EL] = view; view(AL, ELL);
% currAx = axis; axis(currAx);


