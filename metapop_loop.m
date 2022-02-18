function [t,S,I] = metapop_loop(MaxTime,M,n,loci,A,b,f,q,alpha,gamma,muH,muP,rho,mutmatrix,mutlist,G2,Q2,C0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% metapop_loop.m
%
% Carries out the main loop for the metapopulation simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Populations
M1 = find(C0==max(C0),1);
S = zeros(M,n,MaxTime+1);
I = zeros(M,n,n,MaxTime+1);
S(:,1,1) = ones(M,1)*100;
I(M1,1,1,1) = 10;
mutmatrix2 = repmat(permute(mutmatrix,[3,1,2]),[M,1,1]);

for t=1:MaxTime
    SC = S(:,:,t);
    IC = I(:,:,:,t);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Births
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(f>0)
        IH = sum(IC,3);
        N = repmat(sum(SC,2)+sum(IH,2),[1,n,1]);
        births = poiss_rng((A-q*N).*(SC+f*IH));
    else
        N = repmat(sum(SC,2)+sum(sum(IC,3),2),[1,n,1]);
        births = poiss_rng((A-q*N).*(SC));
    end
    births(isnan(births))=0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Host mutations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    B0 = repmat(births,[1,1,n]).*mutmatrix2;
    mutations_out_ind = poiss_rng(muH*B0);
    mutations_out = sum(mutations_out_ind,3);
    list = find(mutations_out>births,1); % Make sure there aren't more mutations than births
    if(~isempty(list))
        list = find(mutations_out>births); % Make sure there aren't more mutations than births
        [list_i,list_j] = ind2sub(size(mutations_out),list); % Make sure there aren't more mutations than births
        for i=1:length(list_i)
            while(mutations_out(list_i(i),list_j(i))>births(list_i(i),list_j(i)))
                r1 = mutlist(list_j(i),ceil(rand*loci));
                mutations_out_ind(list_i(i),list_j(i),r1) = max(0, mutations_out_ind(list_i(i),list_j(i),r1)-1);
                mutations_out(list_i(i),list_j(i)) = sum(mutations_out_ind(list_i(i),list_j(i),:),3);
            end
        end
    end
    mutations_in = permute(sum(mutations_out_ind,2),[1,3,2]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Deaths
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    S_deaths = poiss_rng(b*SC);
    I_deaths = poiss_rng((b+alpha)*IC);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Infections
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    IP = sum(IC,2);
    
    if(isempty(find(IP,1)))
        disp('Pathogen extinct');
        break;
    end

    infections = poiss_rng(inf_function(SC,IP,Q2,mutmatrix,G2,C0,muP,rho,loci,n,M));
    infections_host = sum(infections,3);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Recovery
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(gamma>0)
        recovery = poiss_rng(gamma*IC);
        recovery_host = sum(recovery,3);
    else
        recovery = 0;
        recovery_host = 0;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check that more individuals aren't leaving a class than currently exist
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Susceptible
    S_OUT = S_deaths + infections_host;
    list = find(S_OUT>SC,1);
    if(~isempty(list))
        list = find(S_OUT>SC);
        [list_i,list_j] = ind2sub(size(S_OUT),list); % Make sure there aren't more mutations than births
        for i=1:length(list_i)
            while(S_OUT(list_i(i),list_j(i))>S(list_i(i),list_j(i),t))
                S_OUT(list_i(i),list_j(i)) = S_OUT(list_i(i),list_j(i)) - 1;
                if(rand<S_deaths(list_i(i),list_j(i))/(S_deaths(list_i(i),list_j(i)) + infections_host(list_i(i),list_j(i))))
                    % Update deaths
                    S_deaths(list_i(i),list_j(i)) = S_deaths(list_i(i),list_j(i)) - 1;
                else
                    % Update infections
                    TEMP3 = cumsum(infections(list_i(i),list_j(i),:));                    
                    k1 = find(TEMP3>rand*TEMP3(end),1);
                    infections(list_i(i),list_j(i),k1) = infections(list_i(i),list_j(i),k1) - 1;
                    infections_host(list_i(i),list_j(i)) = infections_host(list_i(i),list_j(i)) - 1;
                end
            end
        end
    end
    
    % Infected
    I_OUT = I_deaths + recovery;
    list = find(I_OUT>IC,1);
    if(~isempty(list))
        list = find(I_OUT>IC);
        [list_i,list_j,list_k] = ind2sub(size(I_OUT),list);
        for i=1:length(list_i)
            while(I_OUT(list_i(i),list_j(i),list_k(i))>I(list_i(i),list_j(i),list_k(i),t))
                I_OUT(list_i(i),list_j(i),list_k(i)) = I_OUT(list_i(i),list_j(i),list_k(i)) - 1;
                if(rand<I_deaths(list_i(i),list_j(i),list_k(i))/(I_deaths(list_i(i),list_j(i),list_k(i)) + recovery(list_i(i),list_j(i),list_k(i))))
                    % Update deaths
                    I_deaths(list_i(i),list_j(i),list_k(i)) = I_deaths(list_i(i),list_j(i),list_k(i)) - 1;
                else
                    % Update recovery
                    recovery(list_i(i),list_j(i),list_k(i)) = recovery(list_i(i),list_j(i),list_k(i)) - 1;
                    recovery_host(list_i(i),list_j(i)) = recovery_host(list_i(i),list_j(i)) - 1;
                end
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update population trackers
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    S(:,:,t+1) = SC + births + mutations_in - mutations_out + recovery_host - S_OUT;
    I(:,:,:,t+1) = IC + infections - I_OUT;
end
t = (0:MaxTime)';