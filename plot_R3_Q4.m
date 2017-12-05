



NumTX = 3; 
N_User = 5;

figure 
index = [0 3 6 9 12];
plot(index,Utility_avg_B(1,:)/K/1e3,'LineWidth',1,'Marker','o','Color',[0 0 0]); hold on;
plot(index,Utility_avg_B(2,:)/K/1e3,'LineWidth',1.5,'Marker','>', 'LineWidth',2,'LineStyle','-.','Color',[0.850980401039124 0.325490206480026 0.0980392172932625]); hold on;
plot(index,Utility_avg_B(3,:)/K/1e3,'LineWidth',1.5,'Marker','v','LineStyle','--','Color',[0 1 0]); hold on;
xlabel('\it{P}^{\rm(b_{0})} \rm[dBm]','FontSize',14,'FontName','Times New Roman');
grid on;
set(gca,'XTick',index)
set(gca,'XTickLabel',{'43','40','37','34','31'},'FontSize',14,'FontName','Times New Roman')
ylabel('Archievable AvgUT [Gbps]','FontSize',14,'FontName','Times New Roman');
h = legend('HomNet [13] - LOS','HomNet [13] - Blockage','HomNet [13] - NLOS');
set(legend,'FontSize',14,'FontName','Times New Roman');

