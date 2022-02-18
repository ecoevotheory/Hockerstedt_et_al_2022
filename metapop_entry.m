function [S,I,poor,well,RANGE,G,X,Y] = metapop_entry(MaxTime,a,b,cH1,cH2,cP1,cP2,q,alpha,beta,gamma,sigma,rho,muH,muP,assortative)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% metapop_entry.m
%
% Entry function for the metapopulation simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check dependencies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(exist('inf_function','file')~=3)
    mex inf_function.c;
end
if(exist('poiss_rng','file')~=3)
    mex poiss_rng.c;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fixed parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = 0;
loci = 4;
n = 2^loci;
gridsize = 6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up metapopulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1=linspace(0,1,gridsize);
[X0,Y0] = meshgrid(x1,x1);
X=X0(:);
Y=Y0(:);
x2=linspace(x1(2),x1(end-1),gridsize+1);
list = find(X>=x2(1) & X<=x2(end) & Y>=x2(1) & Y<=x2(end));
X(list) = [];
Y(list) = [];
[X0,Y0] = meshgrid(x2,x2);
X=[X;X0(:)];
Y=[Y;Y0(:)];
Z = unique([X,Y],'rows');
X=Z(:,1);
Y=Z(:,2);
M = numel(X);
D=sqrt((X*ones(1,M) - ones(M,1)*X').^2 + (Y*ones(1,M) - ones(M,1)*Y').^2);
D(1:(M+1):end)=inf;
G=zeros(M);
G(D<=(1.01/(gridsize-1)))=1;
connections = sum(G,2);

if(assortative==0) % Randomise connections
    
    connections_sum = sum(connections);
    connections_sum_random=0;
    connections_old = connections;
    while(connections_sum_random~=connections_sum)
        G=0*G;
        connections_remaining = sort(connections_old);
        
        for i=1:(M-1)
            count = 0;
            while(connections_remaining(i)>0 && count<1e2)
                connections_cumsum = cumsum(connections_remaining);
                choice = find(rand*connections_cumsum(end)<connections_cumsum,1);
                if(choice~=i && G(i,choice)==0 && G(i,choice)==0)
                    G(i,choice)=1;
                    G(choice,i)=1;
                    connections_remaining(i)=connections_remaining(i)-1;
                    connections_remaining(choice)=connections_remaining(choice)-1;
                end
                count=count+1;
            end
        end
        connections = sum(G,2);
        
        connections_sum_random = sum(connections);
    end
end

% Find poorly and well connected demes
prc = prctile(connections,[20,80]);
poor = (sum(G,2)<=prc(1));
well = (sum(G,2)>=prc(2));

% Temporary matrices for solver
G2 = repmat(G,[1,1,n]);
C0 = repmat(connections,[1,1,n]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Genetics & costs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[RANGE, CH, CP, Q, mutmatrix] = array_values(sigma,loci,n,cH1,cH2,cP1,cP2);
% mutmatrix2 = repmat(permute(mutmatrix,[3,1,2]),[M,1,1]);
mutlist = zeros(n,loci);
Q2 = repmat(permute((beta*repmat(CP',[n,1]).*Q),[3,2,1]),[M,1,1]);
for i=1:n
    mutlist(i,:) = find(mutmatrix(:,i))';
end
A = repmat(a*CH',[M,1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finished setting up, call simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[t,S,I] = metapop_loop(MaxTime,M,n,loci,A,b,f,q,alpha,gamma,muH,muP,rho,mutmatrix,mutlist,G2,Q2,C0);
