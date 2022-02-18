function theory_supp_table_3_analysis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% theory_supp_table_3_analysis.m
% 
% Analyses simulation data for all parameter sets to generate data for
% Supplementary Table 3. Note, requires simulation data in subfolder 'Data'.
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
                    clear disease_prevalence_poor disease_prevalence_well resistance_poor resistance_well infectivity_poor infectivity_well
                    for k=1:length(S_poor(1,1,1,:))
                        count2=count2+1;
                        [disease_prevalence_poor(:,count2),disease_prevalence_well(:,count2),resistance_poor(:,:,count2),resistance_well(:,:,count2),infectivity_poor(:,:,count2),infectivity_well(:,:,count2),~,~] = theory_fig_analysis(double(S_poor(:,:,:,k)),double(S_well(:,:,:,k)),double(IH_poor(:,:,:,k)),double(IH_well(:,:,:,k)),double(IP_poor(:,:,:,k)),double(IP_well(:,:,:,k)),window_full);
                    end
                end
                toc;
                
                count=count+1;
                
                disease_prevalence_well = disease_prevalence_well(window_full,:);
                disease_prevalence_poor = disease_prevalence_poor(window_full,:);
                resistance_well = resistance_well(window_full,:,:);
                resistance_poor = resistance_poor(window_full,:,:);
                infectivity_well = infectivity_well(window_full,:,:);
                infectivity_poor = infectivity_poor(window_full,:,:);
                
                disease_prevalence_well = disease_prevalence_well(:);
                disease_prevalence_poor = disease_prevalence_poor(:);
                resistance_well = resistance_well(:);
                resistance_poor = resistance_poor(:);
                infectivity_well = infectivity_well(:);
                infectivity_poor = infectivity_poor(:);
                
                disease_prevalence_mean(count,1) = mean(disease_prevalence_well,'omitnan');
                disease_prevalence_mean(count,2) = mean(disease_prevalence_poor,'omitnan');                
                resistance_mean(count,1) = mean(resistance_well,'omitnan');
                resistance_mean(count,2) = mean(resistance_poor,'omitnan');
                infectivity_mean(count,1) = mean(infectivity_well,'omitnan');
                infectivity_mean(count,2) = mean(infectivity_poor,'omitnan');
                
                PARAMS(count,:)=[assortative,beta,cH2,cP2];
                
                PROGRESS = count/TOTAL
            end
        end
    end
end
clear PROGRESS TOTAL beta cP2 cH2 count count2 ans i j k i1 j1 k1 assortative R resistance_mean0 resistance_mean_full resistance_mean_trans S_poor S_well I_poor I_well disease_prevalence_poor disease_prevalence_well resistance_poor resistance_well infectivity_poor infectivity_well

% Transmission
for i=1:2
    list = find(PARAMS(:,1)==(i-1));
    NETWORK_D(i,1)=mean(disease_prevalence_mean(list,1)-disease_prevalence_mean(list,2));
    NETWORK_R(i,1)=mean(resistance_mean(list,1)-resistance_mean(list,2));
    NETWORK_I(i,1)=mean(infectivity_mean(list,1)-infectivity_mean(list,2));
    NETWORK_D(i,2)=std(disease_prevalence_mean(list,1)-disease_prevalence_mean(list,2));
    NETWORK_R(i,2)=std(resistance_mean(list,1)-resistance_mean(list,2));
    NETWORK_I(i,2)=std(infectivity_mean(list,1)-infectivity_mean(list,2));
end

% Transmission
for i=1:length(BETA)
    list = find(PARAMS(:,2)==BETA(i));
    BETA_D(i,1)=mean(disease_prevalence_mean(list,1)-disease_prevalence_mean(list,2));
    BETA_R(i,1)=mean(resistance_mean(list,1)-resistance_mean(list,2));
    BETA_I(i,1)=mean(infectivity_mean(list,1)-infectivity_mean(list,2));
    BETA_D(i,2)=std(disease_prevalence_mean(list,1)-disease_prevalence_mean(list,2));
    BETA_R(i,2)=std(resistance_mean(list,1)-resistance_mean(list,2));
    BETA_I(i,2)=std(infectivity_mean(list,1)-infectivity_mean(list,2));
end

% Host costs
for i=1:2
    if(i==1)
        list = find(PARAMS(:,3)<0);
    else
        list = find(PARAMS(:,3)>0);
    end
    HOST_COSTS_D(i,1)=mean(disease_prevalence_mean(list,1)-disease_prevalence_mean(list,2));
    HOST_COSTS_R(i,1)=mean(resistance_mean(list,1)-resistance_mean(list,2));
    HOST_COSTS_I(i,1)=mean(infectivity_mean(list,1)-infectivity_mean(list,2));
    HOST_COSTS_D(i,2)=std(disease_prevalence_mean(list,1)-disease_prevalence_mean(list,2));
    HOST_COSTS_R(i,2)=std(resistance_mean(list,1)-resistance_mean(list,2));
    HOST_COSTS_I(i,2)=std(infectivity_mean(list,1)-infectivity_mean(list,2));
end

% Parasite costs
for i=1:2
    if(i==1)
        list = find(PARAMS(:,4)<0);
    else
        list = find(PARAMS(:,4)>0);
    end
    PAR_COSTS_D(i,1)=mean(disease_prevalence_mean(list,1)-disease_prevalence_mean(list,2));
    PAR_COSTS_R(i,1)=mean(resistance_mean(list,1)-resistance_mean(list,2));
    PAR_COSTS_I(i,1)=mean(infectivity_mean(list,1)-infectivity_mean(list,2));
    PAR_COSTS_D(i,2)=std(disease_prevalence_mean(list,1)-disease_prevalence_mean(list,2));
    PAR_COSTS_R(i,2)=std(resistance_mean(list,1)-resistance_mean(list,2));
    PAR_COSTS_I(i,2)=std(infectivity_mean(list,1)-infectivity_mean(list,2));
end
clear i list

save('theory_supp_table_3_analysis.mat')
