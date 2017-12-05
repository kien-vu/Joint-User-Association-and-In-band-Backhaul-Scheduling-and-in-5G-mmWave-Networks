function [ PL_dB ] = PathLoss( d )
%Path Loss calculattion e
% d : distance in meter
global ntx % In order to determine the pathloss model
%Parameter for LOS at 28 Ghz
alpha3 = 61.4;
beta3 = 2;
sigma3 = 5.8; % dB
% Parameter for NLOS at 28 Ghz
alpha1 = 72;
beta1 = 2.92;
sigma1 = 8.7; % dB

% Parameter for NLOS at 73 Ghz
alpha2 = 82.7;
beta2 = 2.69;
sigma2 = 7.7; % dB
C  = 2500; % MHz carrier frequency
MCL = 70 ; % minimum pathloss for Urban Scenario
% PL_dB =   max(58.83 + 21*log10(C) + 37.6 * log10(d/1000), MCL);
% PL_dB = alpha2 + beta2 * 10 * log10 (d) + sigma2*randn;
%%
PL_dB_LOS = alpha3 + beta3 * 10 * log10 (d) ;%*randn;+ sigma3
PL_dB_NLOS = alpha1 + beta1 * 10 * log10 (d);% + sigma1*randn;

%% Path loss state calculation
beta = 0.006; % Millimeter Wave Cellular Channel Models for System Evaluation

X = [1 0]; % pick state with given probability

if ntx == 1
    PL_dB = PL_dB_LOS; % LOS case
end

if ntx == 2 % Blockage case
    indicator_1 = get_x_from_pmf(X, [ exp(-beta*d)  (1 - exp(-beta*d))]);
    pathloss_dB_Blockage = indicator_1* PL_dB_LOS + (1-indicator_1)* PL_dB_NLOS; % MBS
    PL_dB = pathloss_dB_Blockage;
end

if ntx == 3
    PL_dB = PL_dB_NLOS; % LOS case
end


% % % Parameter for LOS at 73 Ghz
% alpha3 = 69.8;
% beta3 = 2;
% sigma3 = 5.8; % dB
% Parameter for LOS at 2.5 Ghz


% PL_dB = 55.25 + 18.25 * log10 (d) ;% 10 GHZ

%   PL_dB = 55.25 + 18.25 * log10 (d) ;% 10 GHZ

%  PL_dB = 75.85 + 37.3 * log10 (d) ;% Analysis of mm Wave

end
%% Outage state
% a_out = 1/30;
% b_out = 5.2;
% a_los = 1/67.1;
% p_out(k) = max(0, 1 - exp(-a_out*D(k)+b_out));
% p_los(k) = (1-p_out(k))* exp(-a_los*D(k));
% p_nlos(k) = 1 - p_out(k) - p_los(k);
% pathloss_dB_Outage(k) = p_out(k) * 200 + p_los(k)*pathloss_dB_LOS(k) + p_nlos(k)* pathloss_dB_NLOS(k); % MBS
