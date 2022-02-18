function [RANGE, CH, CP, Q, mutmatrix] = array_values(sigma,loci,n,cH1,cH2,cP1,cP2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% array_values.m
%
% Sets values for constant arrays used in the simulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialise genotypes
GENOS = de2bi(0:(n-1));
RANGE = sum(GENOS,2)/loci;

% Calculate infectivity matrix
Q=zeros(n);
for i=1:n
    for j=1:n
        Q(i,j) = sigma^sum((GENOS(i,:)>GENOS(j,:)));
    end
end

% Mutation matrix
mutmatrix = zeros(n);
for i=1:n
    for j=1:n
        if(sum(abs(GENOS(i,:)-GENOS(j,:)))==1)
            mutmatrix(i,j)=1;
        end
    end
end

% Cost functions
if(cH2==0)
    cH2=1e-5;
end
CH = max(0,1-cH1*(1-exp(cH2*RANGE))/(1-exp(cH2)));

if(cP2==0)
    cP2=1e-5;
end
CP = max(0,1-cP1*(1-exp(cP2*RANGE))/(1-exp(cP2)));

