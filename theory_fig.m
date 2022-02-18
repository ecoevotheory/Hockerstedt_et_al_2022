function theory_fig

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% theory_fig.m
% 
% Generates the theory figure for the main text
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
clf
set(gcf,'color','w')
set(gcf,'PaperUnits','centimeters')
xSize = 16; ySize = 9;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 100 xSize*50 ySize*50])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Metapopulation snapshot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;

% Colours for plotting
col1 = [134,83,136]/255;
col2 = [129,179,127]/255;

% Load data
if(~exist('theory_fig_snapshot_data.mat','file'))
    theory_fig_snapshot_data;
end
load('theory_fig_snapshot_data.mat');

% Temporary data
niether = ~well & ~poor;
niether = find(niether);
well = find(well);
poor = find(poor);
RANGETEMP = repmat(permute(RANGE,[2,1]),[length(S(:,1,1)),1,MaxTime+1]);

% Infections by host genotype
IH = permute(sum(I,3),[1,2,4,3]);

% Host frequencies per deme
HTEMP = (S+IH)./(1E-30+repmat(sum(S+IH,2),[1,n,1]));

% Average resistance and infectivity range per deme and connectivity
resistance = permute(sum(HTEMP.*RANGETEMP,2),[3,1,2]);

% Create plot
subplot(2,3,1)
[j,k1,s]=find(G);
plot([X(j) X(k1)]',[Y(j) Y(k1)]','color',0.65*[1,1,1],'linewidth',0.5);
hold on
for k=1:length(niether)
    H1(k)=plot(X(niether(k)),Y(niether(k)),'ok','MarkerFaceColor','w','markersize',5);
end
for k=1:length(well) % well connected
    H2(k)=plot(X(well(k)),Y(well(k)),'ok','MarkerFaceColor',col2,'markersize',5);
end
for k=1:length(poor) % poorly connected
    H3(k)=plot(X(poor(k)),Y(poor(k)),'ok','MarkerFaceColor',col1,'markersize',5);
end
xlim([0,1])
ylim([0,1])
set(gca,'xtick',[],'ytick',[])
axis off

infected = sum(sum(I(:,:,:,end)>0,3),2);
msize = ceil(resistance(end,:)*13)+4;
for k=1:length(H1)
    if(infected(niether(k)))
        set(H1(k),'Marker','s','markersize',msize(niether(k)))
    else
        set(H1(k),'markersize',msize(niether(k)))
    end
end
for k=1:length(H2)
    if(infected(well(k)))
        set(H2(k),'Marker','s','markersize',msize(well(k)))
    else
        set(H2(k),'markersize',msize(well(k)))
    end
end
for k=1:length(H3)
    if(infected(poor(k)))
        set(H3(k),'Marker','s','markersize',msize(poor(k)))
    else
        set(H3(k),'markersize',msize(poor(k)))
    end
end
text(0,1.1,'A','fontsize',16)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cost heatmaps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;

% Load data
if(~exist('theory_fig_heatmaps_analysis.mat','file'))
    theory_fig_heatmaps_analysis;
end
load('theory_fig_heatmaps_analysis.mat');

% Sweep over costs
COSTS_FULL = zeros(length(CH2),length(CP2));
COSTS_TRANS = zeros(length(CH2),length(CP2));
for i=1:length(CH2)
    for j=1:length(CP2)
        list = find(PARAMS(:,3)==CH2(i) & PARAMS(:,4)==CP2(j));
        COSTS_FULL(i,j)=mean(METRIC_FULL(list));
        COSTS_TRANS(i,j)=mean(METRIC_TRANS(list));
    end
end
mn = ceil(100*min(min(COSTS_TRANS(:)),min(COSTS_FULL(:))))/100;
mx = ceil(100*max(max(COSTS_TRANS(:)),max(COSTS_FULL(:))))/100;

subplot(2,3,2)
imagesc(COSTS_TRANS');set(gca,'ydir','normal')
set(gca,'clim',[mn,mx])
text(0.5,4.5*1.1,'B','fontsize',16)

subplot(2,3,3)
imagesc(COSTS_FULL');set(gca,'ydir','normal')
set(gca,'clim',[mn,mx])
text(0.5,4.5*1.1,'C','fontsize',16)

% Define each row of labels.
row1 = {'strong' 'weak' 'weak' 'strong'};
row2 = {'decel.' 'decel.' 'accel.' 'accel.'};

% Combine the rows of labels into a cell array; convert non-strings to strings/character vectors.
% labelArray is an nxm cell array with n-rows of m-tick-lables.
labelArray = [row1; row2];

% To use right or center justification,
labelArray = strjust(pad(labelArray),'center'); % 'left'(default)|'right'|'center

% Combine the rows of labels into individual tick labels
xtickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));

% Assign ticks and labels
% xtickLabels = ['    ',xtickLabels];
for sub=2:3
    subplot(2,3,sub)
    
    ax = gca();
    ax.XTick = 1:4;
    ax.XLim = [0.5,4.5];
    ax.XTickLabel = xtickLabels;
    %     ax.TickLabelInterpreter = 'tex';
    xlabel('Host costs','interpreter','latex','fontsize',16)
    
    ax.YTick = 1:4;
    ax.YLim = [0.5,4.5];
    ax.YTickLabel = xtickLabels;
    ax.YTickLabelRotation = 90;
    ylabel('Parasite costs','interpreter','latex','fontsize',16)
end
map=colormap('bone');
map=flipud(map);
colormap(map);
C=colorbar('eastoutside');
ylabel(C,'Proportion of simulations','interpreter','latex','fontsize',14)
set(C,'ylim',[mn,mx])
drawnow

temp = get(C,'position');
temp_ax_pos = get(gca,'position');
temp(1) = temp(1) + 0.06;
temp(2) = temp_ax_pos(2);
temp(4) = temp_ax_pos(4);
set(C,'position',temp)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dynamics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;

% Colours for plotting
col1 = [134,83,136]/255;
col2 = [129,179,127]/255;

if(~exist('theory_fig_dynamics_data_analysed.mat','file'))
    theory_fig_dynamics_data;
end
load('theory_fig_dynamics_data_analysed.mat')

subplot(2,3,4)
hold on
patch([T' fliplr(T')],[LW_disprev' fliplr(UW_disprev')],col2,'facealpha',0.5,'linestyle','none');
patch([T' fliplr(T')],[LP_disprev' fliplr(UP_disprev')],col1,'facealpha',0.5,'linestyle','none');

plot(T,MW_disprev,'color',col2,'linewidth',2)
plot(T,MP_disprev,'color',col1,'linewidth',2)
xlim([0,T(end)])
ylim([0,0.5])
set(gca,'fontsize',10)
box on
set(gca,'xtick',[0,5000,1e4],'xticklabel',[0,5,10])
xlabel('Time steps $(\times10^4)$','interpreter','latex','fontsize',16)
ylabel('Disease prevalence','interpreter','latex','fontsize',16)
text(0,0.5*1.1,'D','fontsize',16)

subplot(2,3,5)
hold on
patch([T' fliplr(T')],[LW_HR' fliplr(UW_HR')],col2,'facealpha',0.5,'linestyle','none');
patch([T' fliplr(T')],[LP_HR' fliplr(UP_HR')],col1,'facealpha',0.5,'linestyle','none');

plot(T,MW_HR,'color',col2,'linewidth',2)
plot(T,MP_HR,'color',col1,'linewidth',2)
xlim([0,T(end)])
ylim([0,0.6])
set(gca,'fontsize',10)
box on
set(gca,'xtick',[0,5000,1e4],'xticklabel',[0,5,10])
xlabel('Time steps $(\times10^4)$','interpreter','latex','fontsize',16)
ylabel('Resistance','interpreter','latex','fontsize',16)
text(0,0.6*1.1,'E','fontsize',16)

subplot(2,3,6)
hold on
patch([T' fliplr(T')],[LW_PR' fliplr(UW_PR')],col2,'facealpha',0.5,'linestyle','none');
patch([T' fliplr(T')],[LP_PR' fliplr(UP_PR')],col1,'facealpha',0.5,'linestyle','none');

plot(T,MW_PR,'color',col2,'linewidth',2)
plot(T,MP_PR,'color',col1,'linewidth',2)
xlim([0,T(end)])
ylim([0,0.5])
set(gca,'fontsize',10)
box on
set(gca,'xtick',[0,5000,1e4],'xticklabel',[0,5,10])
xlabel('Time steps $(\times10^4)$','interpreter','latex','fontsize',16)
ylabel('Infectivity','interpreter','latex','fontsize',16)
text(0,0.5*1.1,'F','fontsize',16)

% save2pdf('theory_fig.pdf');
