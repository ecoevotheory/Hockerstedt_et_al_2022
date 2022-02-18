function theory_fig_dynamics_data_raw

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% theory_fig_dynamics_data_raw.m
%
% Renames data from the parameter sweep for the main theory figure. 
% Note, requires simulation data in subfolder 'Data'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters
MaxTime = 1e4;
a = 0.11;
b = 0.01;
cH1 = 0.05;
cP1 = 1;
q = (a-b)/100;
loci = 4;
n = 2^loci;
gridsize = 6;
muH = 1e-2;
muP = 1e-3;
alpha = 0.05;
gamma = 0.05;
rho = 5*1e-5;
sigma = 0.2;
simtotal = 200;
T = 0:10:MaxTime;
num_poor = 20;
num_well = 21;
M = 69;
assortative = 1;
beta=0.01;
cH2 = -3;
cP2 = 10;

filename = strcat('Data/sim_',num2str(assortative),'_',num2str(beta),'_',num2str(cH2),'_',num2str(cP2),'.mat');

if(~exist(filename,'file'))
    theory_fig_data;
end
load(filename);

save('theory_fig_dynamics_data_raw.mat')

