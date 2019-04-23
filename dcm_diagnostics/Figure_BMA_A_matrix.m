%load BMA-file first, the specify following three lines:
N_regr= 1; %effect of what regressor you'ld like to view (e.g. 1 = mean connectiivity across all subjects)
N_regions= 8; %specify number of regions
Pp_thresh = 0.90; %specify threshold of connections

%Automatic from here on
figure('units','centimeters','outerposition',[0 0 28 28]);

%load EP and PP
EP=full(vec2mat(BMA.Ep((1+N_regions^2*(N_regr-1)):(N_regions^2*N_regr)),N_regions)');
PP=full(vec2mat(BMA.Pp((1+N_regions^2*(N_regr-1)):(N_regions^2*N_regr)),N_regions)');

%load A-matrix and check what connections correct or wrong
A_matrix=full(vec2mat(BMA.Ep((1+N_regions^2*(N_regr-1)):(N_regions^2*N_regr)),N_regions)');
A_matrix_NaNs=A_matrix;
A_matrix_NaNs(PP<Pp_thresh)=NaN;

%plot matrix
h=imagesc(A_matrix_NaNs);
set(h,'alphadata',~isnan(A_matrix_NaNs))
colorbar;
axis square;

for region=1:length({DCM{1}.xY.name})
    regions(region)={DCM{1}.xY(region).name};
end

if size(A_matrix,1) == size(A_matrix,2) && ...
        size(A_matrix,1) == length({DCM{1}.xY.name})
    set(gca,'YTickLabel',regions,...
        'YTick',1:length({DCM{1}.xY.name}),...
        'XTickLabel',regions,'fontweight','bold','fontsize',14,'XTick',...
        1:length({DCM{1}.xY.name}),'TickLabelInterpreter', 'none');
end
xlabel('\textbf{\underline{From}}','FontSize',22,'Fontweight','bold','Interpreter','latex'); ylabel('\textbf{\underline{To}}','FontSize',22,'Fontweight','bold','Interpreter','latex');
set(gca,'XAxisLocation','top','Tag','connectivity');

%plot title
title_str = ['Connectivity'];
title(title_str,'FontSize',22);

%plot connection strengths
for side1=1:N_regions
    for side2=1:N_regions
        if ~isnan(A_matrix_NaNs(side2,side1))
            text(side1,side2,num2str(round(A_matrix(side2,side1),2)),'HorizontalAlignment','center','Fontweight','bold','Fontsize',18);
        %elseif isnan(A_matrix_NaNs(side2,side1))
        %    text(side1,side2,num2str(round(A_matrix(side2,side1),2)),'HorizontalAlignment','center','Fontweight','normal','Fontsize',10);
        end
    end
end

%%%%%%%%%%%%%%%%%%%
%specify colormap
%%%%%%%%%%%%%%%%%%%
L=length(A_matrix);
indexValue = 0;     % value for which to set a particular color

topColor = [1 0 0];         % color for maximum data value (red = [1 0 0])
indexColor = [1 1 1];       % color for indexed data value (white = [1 1 1])
bottomcolor = [0 0 1];      % color for minimum data value (blue = [0 0 1])

% Calculate where proportionally indexValue lies between minimum and
% maximum values
largest = max(max(A_matrix));
smallest = min(min(A_matrix));
index = L*abs(indexValue-smallest)/(largest-smallest);

% Create color map ranging from bottom color to index color
% Multipling number of points by 100 adds more resolution
customCMap1 = [linspace(bottomcolor(1),indexColor(1),100*index)',...
            linspace(bottomcolor(2),indexColor(2),100*index)',...
            linspace(bottomcolor(3),indexColor(3),100*index)'];

% Create color map ranging from index color to top color
% Multipling number of points by 100 adds more resolution
customCMap2 = [linspace(indexColor(1),topColor(1),100*(L-index))',...
            linspace(indexColor(2),topColor(2),100*(L-index))',...
            linspace(indexColor(3),topColor(3),100*(L-index))'];

customCMap = [customCMap1;customCMap2];  % Combine colormaps

colormap(customCMap)
c=colorbar;
c.Limits=[smallest largest];

set(gcf,'PaperPositionMode','auto');

%direct='';
%Figure_name='';
%saveas(gcf,[direct '/' Figure_name '.bmp']);