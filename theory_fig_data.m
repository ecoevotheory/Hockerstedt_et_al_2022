function theory_fig_data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% theory_fig_data.m
%
% Carries out the main parameter sweep
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

% Variables
BETA = [0.005,0.01];
CH2 = [-10,-3,3,10];
CP2 = [-10,-3,3,10];

% Check data directory exists
if(~exist('Data','dir'))
    mkdir('Data');
end

% Loop over parameters
for assortative=0:1
    for k1=1:length(BETA)
        beta = BETA(k1);
        for j1=1:length(CP2)
            cP2 = CP2(j1);
            for i1=1:length(CH2)
                cH2 = CH2(i1);
                
                filename = strcat('Data/sim_',num2str(assortative),'_',num2str(beta),'_',num2str(cH2),'_',num2str(cP2),'.mat');
                
                if(~exist(filename,'file'))
                    
                    % Set up data arrays - sizes are set for the parameters
                    % used in the main text
                    S_poor = NaN*zeros(num_poor,n,length(T));
                    S_well = NaN*zeros(num_well,n,length(T));
                    IH_poor = NaN*zeros(num_poor,n,length(T));
                    IH_well = NaN*zeros(num_well,n,length(T));
                    IP_poor = NaN*zeros(num_poor,n,length(T));
                    IP_well = NaN*zeros(num_well,n,length(T));
                    poor = NaN*zeros(M,simtotal);
                    well = NaN*zeros(M,simtotal);
                    
                    parfor i=1:simtotal
                        
                        tic;[S,I,poor(:,i),well(:,i),RANGE,~,~,~] = metapop_entry(MaxTime,a,b,cH1,cH2,cP1,cP2,q,alpha,beta,gamma,sigma,rho,muH,muP,assortative);toc;
                        
                        % Add to existing array
                        S_poor(:,:,:,i) = S(poor(:,i)>0,:,T+1);
                        S_well(:,:,:,i) = S(well(:,i)>0,:,T+1);
                        IH_poor(:,:,:,i) = permute(sum(I(poor(:,i)>0,:,:,T+1),3),[1,2,4,3])
                        IH_well(:,:,:,i) = permute(sum(I(well(:,i)>0,:,:,T+1),3),[1,2,4,3])
                        IP_poor(:,:,:,i) = permute(sum(I(poor(:,i)>0,:,:,T+1),2),[1,3,4,2])
                        IP_well(:,:,:,i) = permute(sum(I(well(:,i)>0,:,:,T+1),2),[1,3,4,2])
                        
                    end
                    % Compress data
                    S_poor = uint8(S_poor);
                    S_well = uint8(S_well);
                    IH_poor = uint8(IH_poor);
                    IH_well = uint8(IH_well);
                    IP_poor = uint8(IP_poor);
                    IP_well = uint8(IP_well);
                    poor = poor>0;
                    well = well>0;
                    
                    % Save data
                    save(filename)
                else
                    display('skipping:')
                    filename
                end
            end
        end
    end
end
