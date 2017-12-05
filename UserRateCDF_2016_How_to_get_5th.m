% close all;
% clear CDF_UE1;
% clear CDF_UEB1;
% clear
% clc
% % % Plot the CDF of user rate
% 
% % %% 300 runs
% load MassiveMIMO_OnlyN64_K32Ntx1_N_User12_Bound250_Runs300On08-Nov-2015.mat
% load NewResults_LOS_28GHz_N_N32_K32_NumTX1_N_User_12_Bound2_250_Runs_300On08-Nov-2015.mat


% % %%% 100 runs
% load MassiveMIMO_OnlyN64_K32Ntx1_N_User12_Bound250_Runs100On06-Nov-2015.mat
% load NewResults_LOS_28GHz_N_N32_K32_NumTX1_N_User_12_Bound2_250_Runs_100On06-Nov-2015.mat

% 500 runs
% load NewResults_LOS_28GHz_N_N32_K32_NumTX1_N_User_12_Bound2_250_Runs_500On09-Nov-2015.mat
% load MassiveMIMO_OnlyN64_K32Ntx1_N_User12_Bound250_Runs500On09-Nov-2015.mat
% mean(mean(mean(UserRate-UserRate_A)))

CDF_UE = zeros(K,Iters);
for k=1:K
 for ntx=1:NumTX
    for nu =1:N_User 
        CDF_UE(k,:) = CDF_UE(k,:) + max(UserRate(k,:,ntx,nu),0);
    end
 end
  CDF_UE(k,:) = CDF_UE(k,:)/N_User;
end

CDF_UE_28 = zeros(1,K*Iters);
for k=1:K
  CDF_UE_28(1 + (k-1)*Iters : k*Iters) = CDF_UE(k,:);
end
CDF_UE_28 = sort(CDF_UE_28);
% Get the 5th percentile point
th5per = prctile(CDF_UE_28,5)

 CDF_UE_A = zeros(K,Iters);
 for k=1:K
 for ntx=1:NumTX
    for nu =1:N_User 
        CDF_UE_A(k,:) = CDF_UE_A(k,:) + max(UserRate_A(k,1:Iters,ntx,nu),0);
    end
 end  
  CDF_UE_A(k,:) = CDF_UE_A(k,:)/N_User;

 end

 CDF_UE_B = zeros(K,Iters);
 for k=1:K
 for ntx=1:NumTX
    for nu =1:N_User 
        CDF_UE_B(k,:) = CDF_UE_B(k,:) + max(UserRate_B(k,1:Iters,ntx,nu),0);
    end
 end
 CDF_UE_B(k,:) = CDF_UE_B(k,:)/N_User;
 end
 
CDF_UE_B_28 = zeros(1,K*Iters);
for k=1:K
  CDF_UE_B_28(1 + (k-1)*Iters : k*Iters) = CDF_UE_B(k,:);
end
CDF_UE_B_28 = sort(CDF_UE_B_28);
% Get the 5th percentile point
th5per_B = prctile(CDF_UE_B_28,5)

  CDF_UE_D = zeros(K,Iters);
 for k=1:K
 for ntx=1:NumTX
    for nu =1:N_User 
        CDF_UE_D(k,:) = CDF_UE_D(k,:) + max(UserRate_D(k,1:Iters,ntx,nu),0);
    end
 end
 CDF_UE_D(k,:) = CDF_UE_D(k,:)/N_User;
 end
 
 mean(Utility_avg_higher)/K
 mean(Utility_B_avg_Higher)/K

 
 

% % figure, clf
cdfplot(CDF_UE_28), hold on,
% % cdfplot(LOS_MUE_A/1000), hold on,
cdfplot(CDF_UE_B_28), hold on,
% % cdfplot(LOS_MUE_D/1000), hold on, 
xlabel('UE data rate [Mbps]','FontSize',16,'FontName','Times New Roman');
ylabel('CDF','FontSize',16,'FontName','Times New Roman');
% % title('@ 28GHz-LOS')
% % h = legend('HetNet-Hybrid','Massive MIMO'); %,'HetNet-HD'
% % set(h,'FontSize',14,'FontName','Times New Roman');


% figure, clf
% cdfplot(CDF_UE(8,:)), hold on,
% % cdfplot(LOS_MUE_A/1000), hold on,
% cdfplot(CDF_UE_B(8,:)), hold on,
% % cdfplot(LOS_MUE_D/1000), hold on,