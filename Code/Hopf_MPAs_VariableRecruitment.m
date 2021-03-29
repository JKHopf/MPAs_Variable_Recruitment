% Hopf_MPAs_VariableRecruitment.m
%
% Jess Hopf
% 2021
% email: jess.hopf@gmail.com
%
% This model asks the question of how does the timing of MPA establishment,
% with respect to variability in larval supply, affect:
% 1) the short-term dynamics of fished populations after MPA establishemnt? 
% 2) our ability to detect MPA effects?
%
% Associated scripts & files:
%   - Functions: Functions\Func_ProjPop.m
%   - Plotting: Hopf_Plotting.m
%   - Inital conditions: N1int_0.2Area_v2_20210118.mat
%   - Species parameter values: SppParameterVals.mat
%
% State variables are:
%   - Fish abundance and biomass
%   - Fishery catch (by abundance) and yield (by weight)
%
% Model Assumptions:
%   - periodic variations in larval supply and survival (kelp abundance)
%   - Running as D-indep model
%   - assumes 50:50 sex ratio
%   - Two patch model
%   - Connectivity:
%           1) Sedentary adults & varying scenarios of dispersive larvae
%           2) Moving adults with one dispersive larvae scenario
%
%   - Effort reallocation scenarios (not really using atm):
%       	1) Effort reallocated proportional to the area lost to reserves
%              with varying degrees of increased and decreased effort after
%              reserves are established

%
% Model species:
%   - Kelp Bass
%   - Also has paramters for Blue rockfish (Sebastes mystinus)
%  

%--------------------------------------------------------------------------
clear
addpath('.\Functions')
addpath('.\Data')

% Variable parameters ---------------------------------------------------
	
% Siumlation time:
    % pre-reserves (run for 56 & 19, others are to set u1)
    tF =  56;%  
    % post-reserves
    tR =  19; 

% number of simulations
    nsim = 1000; %
    
% variability for larval recruitment & kelp white noise
    Larv_sd = 0:0.5:1.5;

% shift in recruitment (up to 11 years, beyond that will repeat)
% 0 is bad year, 3 is good larval year when MPA est
    tShift_vec = [0,3];
            
% Species specific parameters:
    % load table
        load SppParameterVals.mat

    % get info for specific species
    SpParas = SppParas(SppParas.sp == 'Kelp bass',:);
    clear SppParas 
    
% fecundity at length paras:
    % Oda et al (1993) for Kelp bass (P. clathratus)
    SpParas.c = 10^-5.57;
    SpParas.d = 2.93;
    
% Juvenile survival values (dictates population growth):
% Chose a value to run
    eig = 1.00; % 0.98; %1.02; %  
    if eig == 1.00
        SpParas.u1 = [3.784243, 3.784243, 3.76, 3.69]; 
    elseif eig == 0.98
        SpParas.u1 = [3.981, 3.985, 3.98, 3.93] ; 
    elseif eig == 1.02
        SpParas.u1 = [3.6, 3.605, 3.594, 3.56]; 
    end          

% number of populations/patches:
    p = 2;
    
% Area = area in reserves (pop1):
    Area = 0.2; %  

% ------------------------------------------------------------------------    
% Pre-allocate matrices:
% N = pop size, Nbio = pop bio, 
% dim = [ageclasses * num pops, sim time, num Larv_sd, num simulations, num loops]
NP = NaN(SpParas.A_max * p, tF+1, length(Larv_sd), nsim, length(tShift_vec));
NPbio = NP; 
NR = NaN(SpParas.A_max * p, tR+1, length(Larv_sd), nsim, length(tShift_vec));
NF = NR; NRbio = NR; NFbio = NR;

% Df = prop die due to fishing 
% (will be same across sims, but may vary withthe looping vec)
% dim = [num Larv_sd, num loops]
DfP = NaN(length(tShift_vec));
DfR = DfP; DfF = DfP;


% ----- LARVAL SUPPLY ----- 
   % individual functions
    % Dims = [sim time, num Larv_sd, num sims]
    R2 = repmat(0.4619.*cos((2*pi/2)*(1:tR+tF+1+max(tShift_vec)))',1,1,nsim);
    R6 = repmat(0.6162.*cos((2*pi/6)*(1:tR+tF+1+max(tShift_vec)))',1,1,nsim);
    Rwn = normrnd(0,repmat(sqrt(0.7018),tR+tF+1+max(tShift_vec),1,nsim)); 
    
    % recruitment function
    RK = (R2 + R6 + Rwn) .* Larv_sd; 
    RKs = exp(RK);
    RKs = RKs./mean(RKs);
 

% ------ Initial conditions --------
    % to get init cons run for a long pre-res time and no var and save output as:
%     N1 = NP(:,end,1,1,1);
%     save 'N1int_0.2Area_v2_age8_20210112.mat'  N1 
 load 'N1int_0.2Area_v2_20210118.mat' % loads N1 to be used in itital conditions
%  N1 = N1.*10000; % testing that the results are population size independent (which they are)


% ---------------------------- Looping --------------------------------
for i = 1:length(tShift_vec)

    % select phase
    tShift = tShift_vec(i);

% ----- PROJECT POPULATIONS -----

 [NP(:,:,:,:,i),NPbio(:,:,:,:,i),DfP(i)] = Func_ProjPop('fished', SpParas,...
                                            repmat(N1,[1,1,length(Larv_sd),nsim]),...
                                            tF, Area, nsim, p,...
                                            RKs, tShift, 1);
                                        
 
% project post-reserves
% with reserves
 [NR(:,:,:,:,i),NRbio(:,:,:,:,i),DfR(i)] = Func_ProjPop('reserve', SpParas,...
                                            NP(:,end,:,:,i), tR, Area, nsim, p,...
                                            RKs(tF+1:end,:,:), tShift, 1);                                           
  
end                                        

% discard the start of the pre-res matrix
NP = NP(:,end-31:end,:,:,:);
NPbio = NPbio(:,end-31:end,:,:,:);

% ages that are fished
AgesF = [zeros(SpParas.Ac-1,1);ones(SpParas.A_max-SpParas.Ac+1,1)];







