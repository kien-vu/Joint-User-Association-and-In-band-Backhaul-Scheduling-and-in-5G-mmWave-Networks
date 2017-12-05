function [] = Systemparameter()
%% -----System Parameters---------
global Axis; % In order change the network size 
Axis = 100;
global Iters; % or time slot
Iters = 300;
global NumTX; % Times to change number of antennas or pathloss model
NumTX = 3; 
global N_User % Times to change number of UEs or others 
N_User = 5;
global ntx
ntx = 1;
global M     % Number of MUE users
M = 8;    
global S     % Number of SCs/SUEs users
S = 4;     
global SCperRow
SCperRow = sqrt(S); % Number of SCs per each row or column
global K    % Total number of Users
K = M + S; 
global N    % Number of MBS Antennas 
N = K;  
global N_s  % Number of Antennas at SCs
N_s = 4;  
global N_au % Number of active users, SC can serve up to N_au users
N_au = 2;  
global F % Number of data streams for spatial multiplexing
F = 1 ; % number of antennas to eliminate the cross-tier interference: Data streams for spartial multiplexing
% Note Update the number of antennas for cross-tier interference after 
% running the user scheduling policy
global P_b0 % Maximum transit power at the macro BSA
P_b0 = 15 + PathLoss(Axis/2) - 84; % 84 in dB at 28GHzglobal Tx   % Maximum transit power at the SC
global Tx
Tx = 38; % dBm
global Tx_Gain
Tx_Gain = 7;% dBi
global Epsilon; % Interference Threshold
Epsilon = 0.001;
global BW; % System bandwidth
BW = 1000; % Ghz
global varphi;
varphi = 2000 * BW; % Define Non-negative weight of penanty function
global Mean
Mean = BW; % The mean of traffic arrival
global scalingFactor
scalingFactor = 1e6;
global ConvP 
ConvP = 0.001; % Convergence tollerence 
% Noise figures: BS 5dB, UE, 9 dB
global C % MHz carrier frequency of conventional cellular network
C  = 2500; 
global epsilon % for the imperfect CSI setting
epsilon = 0.1; % or 0.1 Channel Quality: CSI
global rounding
rounding = 6;
global US_Var

global Bound1 % To ensure the MUEs are located faw way MBS with distance Bound1, and HUEs are located inside the coverage of SCs
global Bound2 % To ensure the SCs are located faw way MBS with distance Bound2
Bound1 = 35 * Axis/500; 
Bound2 = 2.5; %  For dense case, To ensure the MUEs are located faw way SCs with distance Bound2
global Boltzmann_c
Boltzmann_c = 1.381 * 10.^(-23); % (Joule per Kelvin)
global Kelvin_c
Kelvin_c = 290;
global noise_f
noise_f = 9; % dB
noise_f = 10^(noise_f/10);
% noise power = bandwidth × kB × T0 × noise figure (W)% BW in Hz
global noise_p
noise_p = 1000 * BW * Boltzmann_c * Kelvin_c * noise_f;
% noise_p = 10^(noise_p/10)

global coherenceinterval
coherenceinterval = 350; % corresponding to 10 MHz coherence BW, and 35 micro second
global trainingduration % pilot sequence length
trainingduration = 50;

end

%% Some NOTES
% The traffic arrivals follow the Poisson distribution with different mean
% arrival rates. The traffic arrival rates are assumed to be the same for
% all users.
% The maximum traffic arrivals are set as  twice the corressponding mean
% arrival rate. 

% Compute the corresponding channel gain of matrix
% g = abs().^2; if single value uses norm
% Transform the power constraints  using the real function
% The interference threshold depends on the carry frequency, if channel
% model uses the low frequency, then the interference threshold should be
% large enough to avoid remove almost users. In case of high frequence, we
% can make the value of interference threshold to high value.
% Important: Each iteration power is divided equaly for all users. To avoid
% users have to wait for long time to scheduled. queue is quickly stable