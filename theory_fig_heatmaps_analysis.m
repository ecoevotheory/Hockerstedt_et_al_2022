function theory_fig_heatmaps_analysis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% theory_fig_heatmaps_analysis.m
% 
% Analyses simulation data for the theory figure heatmaps (host and
% parasite costs). Note, requires simulation data in subfolder 'Data'.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fixed parameters
MaxTime = 1e4;
a = 0.11;
b = 0.01;
cH1 = 0.05;
cP1 = 1;
q = (a-b)/100;
loci = 4;
n = 2^loci;
gridsize = 6;
M = 69;
muH = 1e-2;
muP = 1e-3;
alpha = 0.05;
gamma = 0.05;
rho = 5*1e-5;
sigma = 0.2;
simtotal = 200;
T = 0:10:MaxTime;

% Variables
BETA = [0.005,0.01];
CH2 = [-10,-3,3,10];
CP2 = [-10,-3,3,10];

% Analysis windows (transient and long-term)
start_full = ceil(length(T)*0.9)+1;
fin_full = length(T);
window_full = start_full:fin_full;
start_trans = ceil(length(T)*0.4)+1;
fin_trans = ceil(length(T)*0.4)+100;
window_trans = start_trans:fin_trans;

% Load data
count=0;
TOTAL = 2*length(BETA)*length(CP2)*length(CH2);
for assortative=0:1
    for k1=1:length(BETA)
        beta = BETA(k1);
        for j1=1:length(CP2)
            cP2 = CP2(j1);
            for i1=1:length(CH2)
                cH2 = CH2(i1);
                
                count2 = 0;
                tic;
                filename = strcat('Data/sim_',num2str(assortative),'_',num2str(beta),'_',num2str(cH2),'_',num2str(cP2),'.mat');
                
                if(exist(filename,'file'))
                    
                    load(filename,'S_poor','S_well','IH_poor','IH_well','IP_poor','IP_well');
                    
                    for k=1:length(S_poor(1,1,1,:))
                        count2=count2+1;
                        [~,~,~,~,~,~,resistance_mean_full(:,count2),~] = theory_fig_analysis(double(S_poor(:,:,:,k)),double(S_well(:,:,:,k)),double(IH_poor(:,:,:,k)),double(IH_well(:,:,:,k)),double(IP_poor(:,:,:,k)),double(IP_well(:,:,:,k)),window_full);
                        [~,~,~,~,~,~,resistance_mean_trans(:,count2),~] = theory_fig_analysis(double(S_poor(:,:,:,k)),double(S_well(:,:,:,k)),double(IH_poor(:,:,:,k)),double(IH_well(:,:,:,k)),double(IP_poor(:,:,:,k)),double(IP_well(:,:,:,k)),window_trans);
                    end
                end
                toc;
                
                count=count+1;
                RM_full(:,:,count) = resistance_mean_full;
                RM_trans(:,:,count) = resistance_mean_trans;
                PARAMS(count,:)=[assortative,beta,cH2,cP2];
                
                resistance_mean0=resistance_mean_full;
                resistance_mean0(:,isnan(sum(resistance_mean_full,1)))=NaN;
                R = resistance_mean0';
                
                METRIC_FULL(count,1) = length(find(R(:,2)>1.05*R(:,1) & R(:,3)>1.05*R(:,4)))/sum(~isnan(R(:,1)));                
                
                resistance_mean0=resistance_mean_trans;
                resistance_mean0(:,isnan(sum(resistance_mean_trans,1)))=NaN;
                R = resistance_mean0';
                
                METRIC_TRANS(count,1) = length(find(R(:,2)>1.05*R(:,1) & R(:,3)>1.05*R(:,4)))/sum(~isnan(R(:,1)));                
                PROGRESS = count/TOTAL
            end
        end
    end
end
clear PROGRESS TOTAL beta cP2 cH2 count count2 ans i j k i1 j1 k1 assortative R resistance_mean0 resistance_mean_full resistance_mean_trans S_poor S_well I_poor I_well

save('theory_fig_heatmaps_analysis.mat')
