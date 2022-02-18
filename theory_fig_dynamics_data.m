function theory_fig_dynamics_data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% theory_fig_dynamics_data.m
% 
% Analyses simulation data for the theory figure dynamics plots.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filename = 'theory_fig_dynamics_data_raw.mat';
if(~exist(filename,'file'))
    theory_fig_dynamics_data_raw;
end
load(filename);

start = ceil(length(S_poor(1,1,:,1))*0.9)+1;
window = start:length(S_poor(1,1,:,1));

count = 0;
for k=1:length(S_poor(1,1,1,:))
    count=count+1;
    tic;[disease_prevalence_poor(:,count),disease_prevalence_well(:,count),resistance_poor(:,:,count),resistance_well(:,:,count),infectivity_poor(:,:,count),infectivity_well(:,:,count),~,~] = theory_fig_analysis(double(S_poor(:,:,:,k)),double(S_well(:,:,:,k)),double(IH_poor(:,:,:,k)),double(IH_well(:,:,:,k)),double(IP_poor(:,:,:,k)),double(IP_well(:,:,:,k)),window);toc;
    PROGRESS = k/length(S_poor(1,1,1,:))
end
clear ans PROGRESS k count S_poor S_well I_poor I_well poor well

MW_disprev = mean(disease_prevalence_well(:,~any(isnan(disease_prevalence_well))),2);
SW_disprev = std(disease_prevalence_well(:,~any(isnan(disease_prevalence_well))),[],2);
MP_disprev = mean(disease_prevalence_poor(:,~any(isnan(disease_prevalence_poor))),2);
SP_disprev = std(disease_prevalence_poor(:,~any(isnan(disease_prevalence_poor))),[],2);

MW_HR = mean(mean(resistance_well,3,'omitnan'),2,'omitnan');
SW_HR = std(mean(resistance_well,2,'omitnan'),[],3,'omitnan');
MP_HR = mean(mean(resistance_poor,3,'omitnan'),2,'omitnan');
SP_HR = std(mean(resistance_poor,2,'omitnan'),[],3,'omitnan');

MW_PR = mean(mean(infectivity_well,3,'omitnan'),2,'omitnan');
SW_PR = std(mean(infectivity_well,2,'omitnan'),[],3,'omitnan');
MP_PR = mean(mean(infectivity_poor,3,'omitnan'),2,'omitnan');
SP_PR = std(mean(infectivity_poor,2,'omitnan'),[],3,'omitnan');

UW_disprev = MW_disprev+SW_disprev;
LW_disprev = max(0,MW_disprev-SW_disprev);
UP_disprev = MP_disprev+SP_disprev;
LP_disprev = max(0,MP_disprev-SP_disprev);

UW_HR = MW_HR+SW_HR;
LW_HR = max(0,MW_HR-SW_HR);
UP_HR = MP_HR+SP_HR;
LP_HR = max(0,MP_HR-SP_HR);

UW_PR = MW_PR+SW_PR;
LW_PR = max(0,MW_PR-SW_PR);
UP_PR = MP_PR+SP_PR;
LP_PR = max(0,MP_PR-SP_PR);

save('theory_fig_dynamics_data_analysed.mat')