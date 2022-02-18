function [disease_prevalence_poor,disease_prevalence_well,resistance_poor,resistance_well,infectivity_poor,infectivity_well,resistance_mean,infectivity_mean] = theory_fig_analysis(S_poor,S_well,IH_poor,IH_well,IP_poor,IP_well,window)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% theory_fig_analysis.m
% 
% Analyses simulation data and returns disease prevalence, resistance and
% infectivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialise genotypes
n = length(S_poor(1,:,1));
loci = log2(n);
GENOS = de2bi(0:(n-1));
RANGE = sum(GENOS,2)/loci;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate disease prevalence per deme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S_POOR_ALL=permute(sum(S_poor,2),[3,1,2]);
S_WELL_ALL=permute(sum(S_well,2),[3,1,2]);
I_POOR_ALL=permute(sum(IH_poor,2),[3,1,2]);
I_WELL_ALL=permute(sum(IH_well,2),[3,1,2]);

disease_prevalence_poor = sum(I_POOR_ALL,2)./sum(S_POOR_ALL+I_POOR_ALL,2); % mean disease prevalence in poorly connected demes
disease_prevalence_well = sum(I_WELL_ALL,2)./sum(S_WELL_ALL+I_WELL_ALL,2); % mean disease prevalence in well connected demes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate average resistance and infectivity range per deme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Temporary matrix
RANGETEMP_POOR = repmat(RANGE',[length(S_poor(:,1,1)),1,length(S_poor(1,1,:))]);
RANGETEMP_WELL = repmat(RANGE',[length(S_well(:,1,1)),1,length(S_well(1,1,:))]);

% Host and parasite frequencies per deme
HTEMP_POOR = (S_poor+IH_poor)./(1E-30+repmat(sum(S_poor+IH_poor,2),[1,n,1]));
HTEMP_WELL = (S_well+IH_well)./(1E-30+repmat(sum(S_well+IH_well,2),[1,n,1]));
PTEMP_POOR = IP_poor./(1E-30+repmat(sum(IP_poor,2),[1,n,1]));
PTEMP_WELL = IP_well./(1E-30+repmat(sum(IP_well,2),[1,n,1]));

% Average resistance and infectivity range per deme and connectivity
resistance_poor = permute(sum(HTEMP_POOR.*RANGETEMP_POOR,2),[3,1,2]);
resistance_well = permute(sum(HTEMP_WELL.*RANGETEMP_WELL,2),[3,1,2]);
infectivity_poor = permute(sum(PTEMP_POOR.*RANGETEMP_POOR,2),[3,1,2]);
infectivity_well = permute(sum(PTEMP_WELL.*RANGETEMP_WELL,2),[3,1,2]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Split data according to connectivity and infection status over
% measurement window
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Temporary matrices
resistance_poor_temp = resistance_poor(window,:);
resistance_well_temp = resistance_well(window,:);
infectivity_poor_temp = infectivity_poor(window,:);
infectivity_well_temp = infectivity_well(window,:);
I_POOR_ALL_TEMP = I_POOR_ALL(window,:);
I_WELL_ALL_TEMP = I_WELL_ALL(window,:);

% Find timepoints where the populations are infected and notinfected
infected_poor = I_POOR_ALL_TEMP>0;
infected_well = I_WELL_ALL_TEMP>0;
notinfected_poor = I_POOR_ALL_TEMP==0;
notinfected_well = I_WELL_ALL_TEMP==0;

% Calculate average resistance and infectivity by connectivity and
% infection status
resistance_mean = [mean(resistance_well_temp(infected_well)),mean(resistance_well_temp(notinfected_well)),mean(resistance_poor_temp(infected_poor)),mean(resistance_poor_temp(notinfected_poor))]';
infectivity_mean = [mean(infectivity_well_temp(infected_well)),mean(infectivity_poor_temp(infected_poor))]';