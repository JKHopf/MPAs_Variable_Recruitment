function [N, Nbio, Df] = Func_ProjPop(scenario, SpParas,...
                                        Ninit, Rtime, Area, nsim, p,...
                                        Rs, tShift, echange)
% Projecting the popualtion over time
% Used for both fully fished or spatially managed pops
%
% Jess Hopf
% 2021
% email: jess.hopf@gmail.com
%
% Inputs:
%   - scenario = 'fished' or 'reserves'
%   - SpParas = single row table of parameters for the species
%   - Ninit = initial population size (matrix size = age-max x nsim)
%   - Rtime = the run time of the projection
%   - Area = area in reserves
%   - nsim = number of simulations
%   - p = number of populations
%   - Rs = matrix of larval peridocity variations over time 
%           (size = at least Rtime x nsim)
%   - tShift = the shift along the Rs
%   - echange = proportional change in effort with reseves 
%               (this is different to fishery squeeze)
%
% Outputs: 
%   - N = population size
%     size is [ageclasses * num pops, sim time, num simulations, num loops]
%   - Nbio = population biomass
%     size is [ageclasses * num pops, sim time, num simulations, num loops]
%   - Df = proportion removed due to fishing
%     size is [num loops, 1]
%
% Associated scripts:
%   - Functions: Func_fecunds.m, Func_Length.m, Func_vecperm.m, 
%                Func_Weights.m
% ------------------------------------------------------------------------
% Pre-allocate matrices:
    N = NaN(SpParas.A_max * p, Rtime, size(Rs,2), nsim);
    Nbio = N;
    
% ------------------------------------------------------------------------    
% Set up demography

% Age, length, weight relationships:
    % Per capita length (cm) at start of age year (L)
    [~,Lengths]= Func_Length((1:SpParas.A_max)',...
                             SpParas.L_inf,SpParas.K,SpParas.A0);

    % Per capita weight at the length at the start of age year (W)
    [~,Weights]=Func_Weights((1:SpParas.A_max)',Lengths,...
                             SpParas.y,SpParas.z);
    
% Fecundity:
    % Per capita eggs produced           
        % Density indep fecundity (per capita) at length at the start 
        % of age year
        % lengths need to be in mm to match co-eff from Dick et al. 2017
        % divide by 2, assuming 50:50 sex ratio
        [~,Fun]=Func_fecunds((1:SpParas.A_max)',Lengths*10,...
                             SpParas.c, SpParas.d);
        Fun = Fun./2;

    % Set ages which reproduce
     Fun(1:(SpParas.A_mat-1))=0;   
           
%      %test
%      Fun = repmat(1/33,33,1);
    
% Fishing mortality:
    switch scenario
        case 'fished'
          mf = SpParas.mf;  
        case 'reserve'
          mf = SpParas.mf/(1-Area) * echange;  
    end
         
    % Proportion caught
    Df = (mf./(SpParas.m+mf))...
                .*(1-exp(-SpParas.m-mf));

% Survival probabilities:
    % survive natural mortality 
    SrU = exp(-SpParas.m);   
    % survive natural & fishing mort 
    % (inc reallocated effort if reserve scenario)
    SrF = exp(-(SpParas.m+mf)); 

% Demographic matrices:  
    % Reefs that are fished 
    RF = diag([repmat(SrU,1,(SpParas.Ac-1)),...
          repmat(SrF,1,SpParas.A_max-SpParas.Ac)],-1);
    RF = repmat(RF,1,1,size(Rs,2),nsim);
    % add fecundity
    RF(1,:,:,:) = repmat(Fun./2,[1,1,size(RF,3:4)]);

    % Reefs that are not fished
    RU = diag(repmat(SrU,1,(SpParas.A_max-1)),-1);
	RU = repmat(RU,1,1,size(Rs,2),nsim);
    % add fecundity
    RU(1,:,:,:) = repmat(Fun./2,[1,1,size(RU,3:4)]);
    
    % Combine pops
    % Dims are [number of populations*max age, number of populations*max age,
    %           num Larv_sd, num sims]
    z = zeros(size(RF));
    switch scenario
        case 'fished'
            B = [RF,z;z,RF];
        case 'reserve'
            B = [RU,z;z,RF]; % [RU,z;z,RU]; %
    end
 
% ------------------------------------------------------------------------
% run over time (inc time-dependent parameters)

% Intial conditions:
N(:,1,:,:) = Ninit;
Nbio(:,1,:,:) = N(:,1,:,:).*[Weights;Weights];


for t = 1:Rtime
    
    % Closed population -----------
    % Add larval variability:
    % Just larva or kelp too
%     varL = Rs(t+tShift+1,:,:).*Ks(t+tShift+1,:,:);
    varL = Rs(t+tShift+1,:,:);
    
    % Larval survival:
    L = exp(-SpParas.u1).*varL;
    L = reshape(L,1,1,[],size(L,3));
    
    % Multiple into projection matrix
    BB = B;
    BB(1,:,:,:) = repmat(BB(1,1:SpParas.A_max,:,:).*L.*Area,1,2);
    BB(SpParas.A_max+1,:,:,:) = repmat(BB(SpParas.A_max+1,SpParas.A_max+1:end,:,:)...
                                          .*L.*(1-Area),1,2);
    
    N(:,t+1,:,:) = pagemtimes(BB,N(:,t,:,:));
    Nbio(:,t+1,:,:) = N(:,t+1,:,:).*[Weights;Weights];  
                           

end


