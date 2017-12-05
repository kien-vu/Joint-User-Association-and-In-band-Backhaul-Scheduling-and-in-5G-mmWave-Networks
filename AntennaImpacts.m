% This file is to plot the impact of number of MBS antennas
% That may show the spectrum efficiency
% Date create: 2016503
% Author: Trung Kien VU

%%
close all;
N = zeros(1,8);
N(1) = 32;
for ntx = 1:8
    if ntx >= 2
%         M = M + 4;
%         S = S + 4;
%         K = M + S;
        if ntx >= 3
            N(ntx) = N(ntx-1) + 64; % start from N
        else
             N(ntx) = N(ntx-1) + 32; % start from N
        end
        if ntx >= 7
             N(ntx) = N(ntx-1) + 100 - 64; % start from N
        end
    end
end
figure, clf
tx = [32    64   128   192   256   320   356   392];
% tx = [1 3 5 7 9 11 13];
BW = 1000 * 41;
 Utility_avg_UA = 1.0e+04 * [ 3.7701    4.2423    4.9781    5.6830    6.2321    6.6131    6.9676    7.4768 ];
 Utility_avg_B = 1.0e+04 * [ 3.8808    3.9573    4.6520    5.1986    5.3292    5.9837    6.3050    7.0206];
% tx = [1 3];

plot(tx, Utility_avg_UA/BW,'b-+','LineWidth',2,...
    'MarkerSize',8), hold on,
% plot(tx, Utility_avg_US/BW,'r-o','LineWidth',2,...
%     'MarkerSize',8), hold on,
plot(tx,Utility_avg_B/BW,'k-s','LineWidth',2,...
    'MarkerSize',8), hold on,
grid;
set(gca,'XTick',tx);
set(gca,'XTickLabel',{'32','64','128','192','256','320','356','392'},'FontSize',16,'FontName','Times New Roman');

xlabel('Number of MBS Antennas, S9, M 32, Ns 32, Axis 250, P BS 43 dBm, varphi 2000','FontSize',16,'FontName','Times New Roman');
ylabel('Spectrum Efficiency [Gbps/Hz]','FontSize',16,'FontName','Times New Roman');
h = legend('HetNet-Hybrid','HomNet'); 
% h = legend('HetNet-Hybrid','HetNet-Hybrid: wo-UA','HomNet'); 
set(h,'FontSize',16,'FontName','Times New Roman');