% Joint Load Balancing and Interference Mitigation in mmWave based Massive MIMO Networks %
% Authors: Trung Kien Vu, 
% Date Update on 2016 March 15
% MATLAB version: 8.4 (R2014b) License host:
%   Username: vkien
%   Host ID: 842b2bc09a54 (eth3)
% CVX Version 2.1, Build 1103 (9714d49)
% Yalmip 

% servers
clc;
clear;
close all;
addpath /home/vkien/yalmip;
addpath /home/vkien/yalmip/extras;
addpath /home/vkien/yalmip/demos;
addpath /home/vkien/yalmip/solvers;
addpath /home/vkien/yalmip/modules;
addpath /home/vkien/yalmip/modules/parametric;
addpath /home/vkien/yalmip/modules/moment;
addpath /home/vkien/yalmip/modules/global;
addpath /home/vkien/yalmip/modules/sos;
addpath /home/vkien/yalmip/operators;
addpath /home/vkien/yalmip/solvers/sedumi;
install_sedumi;
addpath /home/vkien/yalmip/solvers/sdpt3;
install_sdpt3;

 addpath /home/vkien/mosek;
 addpath /home/vkien/mosek/7;
 addpath /home/vkien/mosek/7/toolbox/r2013a;

% 
% % For desktop
clc;
% clear;
% close all;
% addpath K:/yalmip;
% addpath K:/yalmip/extras;
% addpath K:/yalmip/demos;
% addpath K:/yalmip/solvers;
% addpath K:/yalmip/modules;
% addpath K:/yalmip/modules/parametric;
% addpath K:/yalmip/modules/moment;
% addpath K:/yalmip/modules/global;
% addpath K:/yalmip/modules/sos;
% addpath K:/yalmip/operators;
% addpath K:/yalmip/solvers/sedumi;
% install_sedumi;
% addpath K:/yalmip/solvers/sdpt3;
% install_sdpt3;
% 
% addpath 'C:/Program Files/Mosek/7';
% addpath 'C:/Program Files/Mosek/7/toolbox/r2013a';
%% -----System Parameters---------
% Load system parameter
global M  % Number of MUEs
global S  % Number of SCs
global K
global N  % Number of MBS antennas
global F
global epsilon % Channel Quality Index
global BW
global Theta_tilde;
global Epsilon
global sinr_S2u 
global US_Var
global US_Var_A
global scalingFactor
scalingFactor = 1e6;
global random_seed
global ConvP 
% random = rng;% used for generate the same channel parameter for comaprision
global Tx
% Initial Value
global Iters; % or time slot
global NumTX;
global N_s 
global N_au 
global P_b0; % Total transit power at the macro BS
global Axis; % In order change the network size 
global Bound1;

global ntx % In order to determine the pathloss model


%%  Queue Initilization
global MUEnode
global MUEnode_Ad
global HeNBnode
global HUEnode
global THETA
global Theta_S
global ChanG
global Arrival
global Arrival_Ad
global Var
global if_fun
global Tx_Gain
global Convergence
global Loop_max
global N_User
global Mean % the mean of traffic arrival
global varphi; % Define Non-negative weight of penanty function
omega = 1; % weight of proportionally fair utility function
log_index = 0.69314718056;
N0 = -104.5; % dBm for 10 Mhz bandwidth
global noise_p
global coherenceinterval
coherenceinterval = 350; % corresponding to 10 MHz coherence BW, and 35 micro second
global trainingduration % pilot sequence length

global ops; % solver options
ops = sdpsettings('solver','mosek','cachesolvers',1,'verbose',0,'debug',1,'beeponproblem',1,'saveduals',1); % set the interal solver to be mosek
Systemparameter;


%% For Figures
 ConvergencePoints = zeros(Iters,NumTX,N_User);
 NuofMUEs_associatedwithSCs = zeros(Iters,NumTX,N_User);
 Utility_avg_UA = zeros(NumTX,N_User);
 Utility_avg_B = zeros(NumTX,N_User);
 
 Rate_avg = zeros(Iters,NumTX,N_User);
 Rsue_avg = zeros(Iters,NumTX,N_User);
 User_avg = zeros(Iters,NumTX,N_User); % MUE and SUE (min(SC,SUE))
 Q_avg = zeros(Iters,NumTX,N_User);
 Arival_avg = zeros(Iters,NumTX,N_User);
 Y_avg = zeros(Iters,NumTX,N_User);
 C_avg = zeros(Iters,NumTX,N_User);
 D_avg = zeros(Iters,NumTX,N_User);
 Utility =  zeros(Iters,NumTX,N_User);
 
 R_MUE = zeros(Iters,NumTX,N_User);
 R_SC = zeros(Iters,NumTX,N_User);


 Op_Tx_SC = zeros(Iters,NumTX,N_User);


UserS = ones(K+S,Iters,NumTX,N_User); % for users scheduling vector
UA_MUE = ones(S,M,Iters,NumTX,N_User); % for users scheduling vector

Power_Allocate = zeros(K,Iters,NumTX,N_User); % Using CVX
% Data Rate
UserRate = zeros(K,Iters,NumTX,N_User); % Using CVX  % average data rare vector of all users
Rsue = zeros(S,Iters,NumTX,N_User); % The utility function of sue

User_Association_SC = zeros(2*S,Iters,NumTX,N_User); % pre_index for user scheduling
User_Association_MUE = zeros(1+S,M,Iters,NumTX,N_User); % pre_index for user scheduling
SINR_Scheduling_MUE = zeros(1+S,M,Iters,NumTX,N_User); % inside_index for user scheduling
SINR_Scheduling_SC = zeros(2*S,Iters,NumTX,N_User); % inside_index for user scheduling

% After delare the matrixes and variables, lets get this party started
% Initiate the random number generators with a random seed
rng('shuffle'); %  the random number generator based on the current time so that rand, randi, and randn produce a different sequence of numbers after each time you call rng.

random_seed = rng; %  seeds the random number generator using the nonnegative integer seed so that rand, randi, and randn produce a predictable sequence of numbers.

Arrival = poissrnd(Mean,K,Iters)'; 
Arrival_Ad = ones(Iters,S);
% Initiate the random number generators with a random seed

fprintf('\n Reviewer 3 Question 3\n');

rng(random_seed);
US_Var = ones(1,K+S);
[UserRate_B, Utility_B, Utility_avg_B, Power_Allocate_B ] = MassiveMIMO_Only();
Utility_avg_B