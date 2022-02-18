function theory_fig_snapshot_data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fig_snapshot_data.m
%
% Generates simulation data for the metapopulation snapshot figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reset random number generator
rng default

% Parameters
MaxTime = 1e4;
a = 0.11;
b = 0.01;
cH1 = 0.05;
cP1 = 1;
cH2 = -3;
cP2 = 10;
q = (a-b)/100;
muH = 1e-2;
muP = 1e-3;
alpha = 0.05;
beta = 0.01;
gamma = 0.05;
rho = 5*1e-5;
sigma = 0.2;
assortative = 1;   

% Carry out simulation
tic;[S,I,poor,well,RANGE,G,X,Y] = metapop_entry(MaxTime,a,b,cH1,cH2,cP1,cP2,q,alpha,beta,gamma,sigma,rho,muH,muP,assortative);toc;

save('theory_fig_snapshot_data.mat');