function fitness_cost_fig

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fitness_cost_fig.m
% 
% Generates the fitness cost figure for the supplementary material
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up fixed parameters
L = 4;
c1 = 1;

% Set up variables
z = linspace(0,L,101);
c2 = [-10,-3,3,10];
[Z,C2] = meshgrid(z,c2);

% Calculate fitness costs
cost = c1*(1-exp(C2.*Z/L))./(1-exp(C2));

% Set up figure parameters
lab_locs = [0.7,1.2,2.7,3.2];
spacing = [0.1,0.1,-0.8,-0.8];
lab = '$c_k^2=';

% Create figure
figure(1)
clf
set(gcf,'color','w')
set(gcf,'PaperUnits','centimeters')
xSize = 5; ySize = 5;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 100 xSize*50 ySize*50])

% Update line styles
h=plot(z,cost,'color','k','linewidth',2);
set(h(1),'color','r')
set(h(2),'color','r','linestyle','--')
set(h(3),'linestyle','--')
set(h(3),'linestyle','--')

% Add text labels
for i=1:length(c2)
    loc=find(z>lab_locs(i),1);
    text(z(loc)+spacing(i),cost(i,loc),strcat(lab,num2str(c2(i)),'$'),'interpreter','latex');
end

% Update axes labels
xlabel('Number of resistance or infectivity alleles','interpreter','latex');
set(gca,'ytick',0);
text(-0.3,0.2,'$\frac{c_k^1}{5}$','interpreter','latex','fontsize',12);
text(-0.3,0.4,'$\frac{2c_k^1}{5}$','interpreter','latex','fontsize',12);
text(-0.3,0.6,'$\frac{3c_k^1}{5}$','interpreter','latex','fontsize',12);
text(-0.3,0.8,'$\frac{4c_k^1}{5}$','interpreter','latex','fontsize',12);
text(-0.3,1,'$c_k^1$','interpreter','latex','fontsize',12);
drawnow
y1=ylabel('Fitness cost, $c_H(i)$ or $c_P(j)$','interpreter','latex');
temp = get(y1,'position');
temp(1) = temp(1)-0.1;
set(y1,'position',temp);

% save2pdf('fitness_cost_fig.pdf');