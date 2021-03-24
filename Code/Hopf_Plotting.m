% Hopf_Plotting.m
%
% Jess Hopf
% 2021
% email: jess.hopf@gmail.com
%
% This script creates various plots of model outputs from 
% Hopf_MPAs_VariableRecruitment.m
%

%% simulation lines & mean abundance (normalised)

NRt = [NP(:,1:end-1,:,:,:),NR];

% Larv_sd: [0,0.5,1,1.5]
lsp = 3;
% nsim sample
sims =  randsample(1000,100);
% min age sampled (fished age (Ac) = 8, age at mat = 4 (A_mat))
Samp_age = 2;% SpParas.Ac; %     

tplot = size(NP,2)-6:size(NP,2)+10;
xp = repmat((-6:10)',1,size(sims,2));

figure 
hold on
vc = inferno(6);
vc_vec = [vc(5,:);vc(2,:)];

for i = 1:2
subplot(3,1,1)
hold on
plot(xp,...
    squeeze(NRt(1,tplot,lsp,sims,i)./mean(NRt(1,:,lsp,sims,i),2)),...
    'Color', [vc_vec(i,:),0.1]);
errorbar(xp,...
    squeeze(mean(NRt(1,tplot,lsp,:,i)./mean(NRt(1,:,lsp,:,i),2),4)),...
    squeeze(std(NRt(1,tplot,lsp,:,i)./mean(NRt(1,:,lsp,:,i),2),[],4)),...
    'Color', vc_vec(i,:).*0.5,'LineWidth',1.5)
plot(xp,...
    squeeze(NRt(1,tplot,1,1,i)./mean(NRt(1,:,1,1,i),2)),...
    'Color', [0.7,0.7,0.7],'LineWidth',1.5);
xline(0,'--k')
ylabel("Recrtuiment")
title("(pop growth = " + eig + "; Larv sd = " + Larv_sd(lsp) + "; Age = " + Samp_age +  " )")
ylim([0,10])

subplot(3,1,2)
hold on
plot(xp,...
    squeeze(sum(NRt(Samp_age:33,tplot,lsp,sims,i),1)./sum(NRt(Samp_age:33,size(NP,2),lsp,sims,i),1)),...
    'Color', [vc_vec(i,:),0.1])
errorbar(xp,...
    squeeze(mean(sum(NRt(Samp_age:33,tplot,lsp,:,i),1)./sum(NRt(Samp_age:33,size(NP,2),lsp,:,i),1),4)),...
    squeeze(std(sum(NRt(Samp_age:33,tplot,lsp,:,i),1)./sum(NRt(Samp_age:33,size(NP,2),lsp,:,i),1),[],4)),...
    'Color', vc_vec(i,:).*0.5,'LineWidth',1.5)
plot(xp,...
    squeeze(sum(NRt(Samp_age:33,tplot,1,1,i),1)./sum(NRt(Samp_age:33,size(NP,2),1,1,i),1)),...
    'r','LineWidth',1.5)
xline(0,'--k')
yline(1,':k')
ylabel("Reserve Abund")
ylim([0,5])

subplot(3,1,3)
hold on
plot(xp,...
    squeeze(sum(NRt((33+Samp_age):66,tplot,lsp,sims,i),1)./sum(NRt((33+Samp_age):66,size(NP,2),lsp,sims,i),1)),...
    'Color', [vc_vec(i,:),0.1])
errorbar(xp,...
    squeeze(mean(sum(NRt((33+Samp_age):66,tplot,lsp,:,i),1)./sum(NRt((33+Samp_age):66,size(NP,2),lsp,:,i),1),4)),...
    squeeze(std(sum(NRt((33+Samp_age):66,tplot,lsp,:,i),1)./sum(NRt((33+Samp_age):66,size(NP,2),lsp,:,i),1),[],4)),...
    'Color', vc_vec(i,:).*0.5,'LineWidth',1.5)
plot(xp,...
    squeeze(sum(NRt((33+Samp_age):66,tplot,1,1,i),1)./sum(NRt((33+Samp_age):66,size(NP,2),1,1,i),1)),...
   'r','LineWidth',1.5)
xline(0,'--k')
yline(1,':k')
ylabel("Fished Abund")
xlabel("Years after MPA")
ylim([0,5])
end

%% AUC 
clear AUC1 AUC2 AUC3

NRt = [NP(:,1:end-1,:,:,:),NR];

samplet = size(NP,2)+(1:10);
larv = 2:4; 
Samp_age =  2;%SpParas.Ac; %  

    FishAge = zeros(SpParas.A_max*2,1);
    ResAge = FishAge;
    FishAge(SpParas.A_max+Samp_age:SpParas.A_max*2) = 1./(1-Area);
    ResAge(Samp_age:SpParas.A_max) = 1./Area;

figure
hold on
vc = inferno(6);
colororder ([vc(5,:);vc(2,:)])
lsty = ["--","-","-.",":"];

for h = larv 
for i = 1:2 % [1,4] 
for j = 1:length(samplet)
    
    ts = samplet(j);
    ts1 = size(NP,2);
    
    % 3-yr means
%     ts = ts-1:ts+1; 
%     ts1 = ts1-1:ts1+1;
    
    % ROC hypo: X1 (fished pop) < X2 (res pop)
    [FP,TP] = ROC(mean(sum(NRt(:,ts,h,:,i).*FishAge),2),...
                  mean(sum(NRt(:,ts,h,:,i).*ResAge),2));
    [~,Tind] = unique(FP);
    AUC1(i,j) = -trapz(FP,TP);
    
    [FP,TP] = ROC(mean(sum(NRt(:,ts1,h,:,i).*ResAge),2),...
                  mean(sum(NRt(:,ts,h,:,i).*ResAge),2)); 
    [~,Tind] = unique(FP);
    AUC2(i,j) = -trapz(FP,TP); 
    
    [FP,TP] = ROC(mean(sum(NRt(:,ts,h,:,i).*FishAge),2)...
                    ./mean(sum(NRt(:,ts1,h,:,i).*FishAge),2),...
                  mean(sum(NRt(:,ts,h,:,i).*ResAge),2)...
                    ./mean(sum(NRt(:,ts1,h,:,i).*ResAge),2));
    [~,Tind] = unique(FP);
    AUC3(i,j) = -trapz(FP,TP); 

% ROC plot:    
%   - as the plot bows towards the upper left, it becomes easier to detect
%   a difference (the distributions decrease overlap)
%   - a straight diagonal line means there is a lot of overlap between the
%   distributions, and it's impossible to detect a difference
% plot(FP,TP)
% xlabel('False Positive')
% ylabel('True Positive')
% plot([0,1],[0,1],":k",'linewidth',1)

end
end
subplot(1,3,1)
hold on
plot(samplet-size(NP,2),AUC1','linewidth',1,'linestyle',lsty(find(larv == h)))
title('Outside < Inside')
ylim([0.4,1])
xlim([0 11])
yline(0.5, ":k")

subplot(1,3,2)
hold on
plot(samplet-size(NP,2),AUC2','linewidth',1,'linestyle',lsty(find(larv == h)))
title('Before (in) < After (in)')
ylim([0.4,1])
xlim([0 11])
yline(0.5, ":k")
ylabel('AUC')

subplot(1,3,3)
hold on
plot(samplet-size(NP,2),AUC3','linewidth',1,'linestyle',lsty(find(larv == h)))
title('Change Outside < Change Inside')
ylim([0.4,1])
xlim([0 11])
yline(0.5, ":k")
xlabel('Sample (yrs after est)')
end

lgd = legend("Larv sd: "+string(repelem(Larv_sd(larv),2))+", Phase: "+string([1,4,1,4,1,4]-1));
% lgd = legend(string([0,3]));
lgd.Title.String = "Phase (pop growth = " + eig + "; Age = " + Samp_age +  " )";

%% AUC for different min sample age
clear AUC1 AUC2 AUC3
samplet = size(NP,2)+(1:10); 
larv = 3;
phase = 1;
Samp_age_range = 1:8;

NRt = [NP(:,1:end-1,:,:,:),NR]; 

figure
hold on
vc = plasma(9);
colororder(flip(vc(1:8,:)))
lsty = ["-","--","-.",":"];

for h = larv 
for i = 1:length(Samp_age_range) 
    FishAge = zeros(SpParas.A_max*2,1);
    ResAge = FishAge;
    FishAge(SpParas.A_max+Samp_age_range(i):SpParas.A_max*2) = 1./(1-Area);
    ResAge(Samp_age_range(i):SpParas.A_max) = 1./Area;
    
for j = 1:length(samplet)  
    % ROC hypo: X1 (fished pop) < X2 (res pop)
    [FP,TP] = ROC(sum(NRt(:,samplet(j),h,:,phase).*FishAge) ,...
                  sum(NRt(:,samplet(j),h,:,phase).*ResAge));
    [~,Tind] = unique(FP);
    AUC1(i,j) = -trapz(FP,TP); % mean(TP(Tind));
    
    [FP,TP] = ROC(sum(NRt(:,size(NP,2),h,:,phase).*ResAge) ,...
                  sum(NRt(:,samplet(j),h,:,phase).*ResAge)); 
    [~,Tind] = unique(FP);
    AUC2(i,j) = -trapz(FP,TP); % mean(TP(Tind));
    
    [FP,TP] = ROC(sum(NRt(:,samplet(j),h,:,phase).*FishAge)...
                    ./sum(NRt(:,size(NP,2),h,:,phase).*FishAge),...
                  sum(NRt(:,samplet(j),h,:,phase).*ResAge)...
                    ./sum(NRt(:,size(NP,2),h,:,phase).*ResAge));
    [~,Tind] = unique(FP);
    AUC3(i,j) = -trapz(FP,TP); % mean(TP(Tind));
    
%     % 3-yr means
%     [FP,TP] = ROC(mean(sum(NRt(:,samplet(j)-1:samplet(j)+1,h,:,phase).*FishAge)) ,...
%                   mean(sum(NRt(:,samplet(j)-1:samplet(j)+1,h,:,phase).*ResAge)));
%     [~,Tind] = unique(FP);              
% 	AUC1(i,j) = -trapz(FP,TP); % mean(TP(Tind));
%                   
%     [FP,TP] = ROC(mean(sum(NRt(:,size(NP,2)-1:size(NP,2)+1,h,:,phase).*ResAge)),...
%                   mean(sum(NRt(:,samplet(j)-1:samplet(j)+1,h,:,phase).*ResAge))); 
% 	[~,Tind] = unique(FP);
%     AUC2(i,j) = -trapz(FP,TP); % mean(TP(Tind));
%     
%     [FP,TP] = ROC(mean(sum(NRt(:,samplet(j)-1:samplet(j)+1,h,:,phase).*FishAge))...
%                     ./mean(sum(NRt(:,size(NP,2)-1:size(NP,2)+1,h,:,phase).*FishAge)),...
%                   mean(sum(NRt(:,samplet(j)-1:samplet(j)+1,h,:,phase).*ResAge))...
%                     ./mean(sum(NRt(:,size(NP,2)-1:size(NP,2)+1,h,:,phase).*ResAge)));
%     [~,Tind] = unique(FP);
%     AUC3(i,j) = -trapz(FP,TP); % mean(TP(Tind));                

end
end
subplot(1,3,1)
hold on
plot(samplet-size(NP,2),AUC1','linewidth',1,'linestyle',lsty(find(larv == h)))
title('Outside < Inside')
ylim([0.4,1])
xlim([0 11])
yline(0.5, ":k")

subplot(1,3,2)
hold on
plot(samplet-size(NP,2),AUC2','linewidth',1,'linestyle',lsty(find(larv == h)))
title('Before (in) < After (in)')
ylim([0.4,1])
xlim([0 11])
yline(0.5, ":k")
ylabel('AUC')

subplot(1,3,3)
hold on
plot(samplet-size(NP,2),AUC3','linewidth',1,'linestyle',lsty(find(larv == h)))
title('Change Outside < Change Inside')
ylim([0.4,1])
xlim([0 11])
yline(0.5, ":k")
xlabel('Sample (yrs after est)')
end

[~,Lengths] = Func_Length((1:SpParas.A_max)',...
                             SpParas.L_inf,SpParas.K,SpParas.A0);

lgd = legend(string(Samp_age_range) + " yrs; ~" + string(round(Lengths(Samp_age_range)))' + " cm");
lgd.Title.String = "Min Age Sampled (pop growth = " + eig + "; Phase = " + tShift_vec(phase) +  " )";

%% Total Pop: simulation lines & mean abundance (normalised)

NRt = [NP(:,1:end-1,:,:,:),NR]; % 

% larval_sd: [0,0.5,1,1.5]
lsp = 3;
% nsim sample
sims = randsample(1000,100);

tplot = size(NP,2)-6:size(NP,2)+10; 
xp = repmat((-6:10)',1,size(sims,2));

% sampling ages
ages1 = [2:33,(34+1):66];
ages2 = [SpParas.Ac:33,(34+SpParas.Ac):66];

figure
hold on
vc = inferno(6);
vc_vec = [vc(5,:);vc(2,:)];

for i = 1:2
subplot(3,1,1)
hold on
plot(xp,...
    squeeze(NRt(34,tplot,lsp,sims,i)./mean(NRt(34,:,lsp,sims,i),2)),...
    'Color', [vc_vec(i,:),0.1]);
errorbar(xp,...
    squeeze(mean(NRt(34,tplot,lsp,:,i)./mean(NRt(34,:,lsp,:,i),2),4)),...
    squeeze(std(NRt(34,tplot,lsp,:,i)./mean(NRt(34,:,lsp,:,i),2),[],4)),...
    'Color', vc_vec(i,:).*0.5,'LineWidth',1.5)
plot(xp,...
    squeeze(NRt(34,tplot,1,1,i)./mean(NRt(34,:,1,1,i),2)),...
    'Color', [0.7,0.7,0.7],'LineWidth',1.5);
xline(0,'--k')
ylabel("Recrtuiment")
title("(pop growth = " + eig + "; Larv sd = " + Larv_sd(lsp) +  " )")
ylim([0,10])

subplot(3,1,2)
hold on
plot(xp,...
    squeeze(sum(NRt(ages1,tplot,lsp,sims,i),1)./sum(NRt(ages1,size(NP,2),lsp,sims,i),1)),...
    'Color', [vc_vec(i,:),0.1])
errorbar(xp,...
    squeeze(mean(sum(NRt(ages1,tplot,lsp,:,i),1)./sum(NRt(ages1,size(NP,2),lsp,:,i),1),4)),...
    squeeze(std(sum(NRt(ages1,tplot,lsp,:,i),1)./sum(NRt(ages1,size(NP,2),lsp,:,i),1),[],4)),...
    'Color', vc_vec(i,:).*0.5,'LineWidth',1.5)
plot(xp,...
    squeeze(sum(NRt(ages1,tplot,1,1,i),1)./sum(NRt(ages1,size(NP,2),1,1,i),1)),...
    'r','LineWidth',1.5)
xline(0,'--k')
yline(1,':k')
ylabel("Total Abund: 2+")
ylim([0,5])

subplot(3,1,3)
hold on
plot(xp,...
    squeeze(sum(NRt(ages2,tplot,lsp,sims,i),1)./sum(NRt(ages2,size(NP,2),lsp,sims,i),1)),...
    'Color', [vc_vec(i,:),0.1])
errorbar(xp,...
    squeeze(mean(sum(NRt(ages2,tplot,lsp,:,i),1)./sum(NRt(ages2,size(NP,2),lsp,:,i),1),4)),...
    squeeze(std(sum(NRt(ages2,tplot,lsp,:,i),1)./sum(NRt(ages2,size(NP,2),lsp,:,i),1),[],4)),...
    'Color', vc_vec(i,:).*0.5,'LineWidth',1.5)
plot(xp,...
    squeeze(sum(NRt(ages2,tplot,1,1,i),1)./sum(NRt(ages2,size(NP,2),1,1,i),1)),...
    'r','LineWidth',1.5)
xline(0,'--k')
yline(1,':k')
ylabel("Total Abund: 6+")
ylim([0,5])

end

%% Total Pop: AUC

clear AUC1 AUC2 AUC3

NRt = [NP(:,1:end-1,:,:,:),NR];

samplet = size(NP,2)+(1:10);
larv = 2:4;
Samp_age = 2;% SpParas.Ac; % 

    FishAge = zeros(SpParas.A_max*2,1);
    ResAge = FishAge;
    FishAge(SpParas.A_max+Samp_age:SpParas.A_max*2) = 1./(1-Area);
    ResAge(Samp_age:SpParas.A_max) = 1./Area;

figure
hold on
vc = inferno(6);
colororder ([vc(5,:);vc(2,:)])
lsty = ["--","-","-.",":"];

for h = larv 
for i = 1:2 % [1,4] 
for j = 1:length(samplet)
    % ROC hypo: X1 (after) < X2 (before)
    [FP,TP] = ROC(sum(NRt(:,size(NP,2),h,:,i).*FishAge) + sum(NRt(:,size(NP,2),h,:,i).*ResAge) ,...
                  sum(NRt(:,samplet(j),h,:,i).*FishAge) + sum(NRt(:,samplet(j),h,:,i).*ResAge));
%     [FP,TP] = ROC(sum(NRt(:,samplet(j),h,:,i).*FishAge) + sum(NRt(:,samplet(j),h,:,i).*ResAge),...
%                   sum(NRt(:,size(NP,2),h,:,i).*FishAge) + sum(NRt(:,size(NP,2),h,:,i).*ResAge));

    [~,Tind] = unique(FP);
    AUC1(i,j) = -trapz(FP,TP); %mean(TP(Tind));
end
end
plot(samplet-size(NP,2),AUC1','linewidth',1,'linestyle',lsty(find(larv == h)))
title('Before < After')
ylim([0,1])
xlim([0 11])
yline(0.5, ":k")
ylabel('AUC')
xlabel('Sample (yrs after est)')
end

lgd = legend("Larv sd: "+string(repelem(Larv_sd(larv),2))+", Phase: "+string([1,4,1,4,1,4]-1));
lgd.Title.String = "Phase (pop growth = " + eig + "; Age = " + Samp_age +  " )";




