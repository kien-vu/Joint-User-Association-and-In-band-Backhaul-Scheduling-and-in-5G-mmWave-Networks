
load FigurePlotsN6_K10Runs30Lya4000.mat;
% load MassiveMIMO_OnlyN32_K36_28-May-2015.mat;
% load MassiveMIMO_Wired_N32_K36_28-May-2015.mat;
% load MassiveMIMO_WirelessN32_K36_28-May-2015.mat;
 Iters = 30;
 NumTX =1;
%% Draw figures
% Iters = 20
figure(1), clf
axis([0 Iters 0 5])
% for ntx=1:NumTX
plot( 1:Iters,User_avg_B(:,1), 'k-+'), hold on,
plot( 1:Iters,User_avg(:,1), 'b-*'), hold on,
plot( 1:Iters,User_avg_A(:,1), 'r-.'), hold on,
% end
xlabel('Iters)','FontSize',14,'FontName','Times New Roman');
ylabel('Average Data Rate [Mbps]','FontSize',14,'FontName','Times New Roman');
legend('Massive MIMO','HetNet-Wireless','HetNet-Wired');
set(legend,'FontSize',14,'FontName','Times New Roman');


figure(3), clf
axis([0 Iters 0 5])
% for ntx=1:NumTX
plot( 1:Iters,R_MUE(:,1), 'b-*'), hold on,
plot( 1:Iters,User_avg_A(:,1), 'r-.'), hold on,
% end
xlabel('Iters','FontSize',14,'FontName','Times New Roman');
ylabel('Average MUE Data Rate [Mbps]','FontSize',14,'FontName','Times New Roman');
legend('HetNet-Wireless','HetNet-Wired');
set(legend,'FontSize',14,'FontName','Times New Roman');

figure(4), clf
axis([0 Iters 0 5])
% for ntx=1:1
plot( [1:Iters],Rsue_avg(:,1), 'b-*', [1:Iters],Rsue_avg_A(:,1), 'k-h'), hold on,
% end
xlabel('Iters','FontSize',14,'FontName','Times New Roman');
ylabel('Average SUE Data Rate [Mbps]','FontSize',14,'FontName','Times New Roman');
legend('HetNet-Wireless','HetNet-Wired');
set(legend,'FontSize',14,'FontName','Times New Roman');
% 


figure, clf
plot( [1:Iters],Utility_B(:,1),'k-h', 'LineWidth',1 ), hold on,
plot( [1:Iters],Utility(:,1), 'b-*', 'LineWidth',1 ), hold on,
plot( [1:Iters],Utility_A(:,1), 'r-.', 'LineWidth',1 ), hold on,
xlabel('Iters','FontSize',14,'FontName','Times New Roman');
ylabel('Total Network Utility [Mbps]','FontSize',14,'FontName','Times New Roman');
% legend('Massive MIMO','HetNet-Wireless');
legend('Massive MIMO','HetNet-Wireless','HetNet-Wired');
set(legend,'FontSize',14,'FontName','Times New Roman');