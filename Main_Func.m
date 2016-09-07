% Joint Load Balancing and Interference Mitigation in 5G HetNets & Joint
% In-Band Scheduling and Interference in 5G Heterogeneous Network
% Authors: Trung Kien Vu,  
% Date Update on 2016 September 06
% MATLAB version: 8.4 (R2014b) License host:
%   Username: vkien
%   Host ID: 842b2bc09a54 (eth3)
% CVX Version 2.1, Build 1103 (9714d49)
% YALMIP R20150918 (Release announcement)
% Do not forget that you typically need solvers besides those installed by default in MATLAB

% % Add solvers manually if you run on your laptop
% addpath /home/vkien/yalmip;
% addpath /home/vkien/yalmip/extras;
% addpath /home/vkien/yalmip/demos;
% addpath /home/vkien/yalmip/solvers;
% addpath /home/vkien/yalmip/modules;
% addpath /home/vkien/yalmip/modules/parametric;
% addpath /home/vkien/yalmip/modules/moment;
% addpath /home/vkien/yalmip/modules/global;
% addpath /home/vkien/yalmip/modules/sos;
% addpath /home/vkien/yalmip/operators;
% addpath /home/vkien/yalmip/solvers/sedumi;
% addpath /home/vkien/yalmip/SeDuMi;
% install_sedumi;
% addpath /home/vkien/yalmip/sdpt3;
% install_sdpt3;

% 
% % For desktop
clc;
clear;
close all;
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
% addpath K:/yalmip/solvers/mosek;

%% -----System Parameters---------
% Load system parameter
global M  % Number of MUEs
global S  % Number of SCs
global K  % Total Number of Users
global N  % Number of MBS antennas
global F  % Number of Stream Data for each User
global epsilon % Channel Quality Index
global BW  % System Bandwidth
global Theta_tilde; % the spatial channel correlation matrix
global scalingFactor
scalingFactor = 1e6;
global random_seed
global Iters; % or time slot

%%  Location, Queue and Traffic Initilization
global MUEnode
global MUEnode_Ad
global HeNBnode
global HUEnode
global THETA  % the spatial channel correlation matrix MUEs
global Theta_S % the spatial channel correlation matrix SCs
global ChanG % Channel generization
global Arrival % Arrival Rate
global Arrival_Ad
global Mean % the mean of traffic arrival
global varphi; % Define Non-negative weight of penanty function
global coherenceinterval
coherenceinterval = 350; % corresponding to 10 MHz coherence BW, and 35 micro second
global trainingduration % pilot sequence length

global ops; % solver options for Yalmip
% ops = sdpsettings('solver','sedumi','verbose',0,'sedumi.eps',1e-5);
% ops = sdpsettings('solver','sdpt3','verbose',0); % set the interal solver to be SDPT3
% ops = sdpsettings('solver','sdpt3','sdpt3.maxit',300,'cachesolvers',1,'verbose',0); % set the interal solver to be SDPT3
% ops = sdpsettings('solver','sedumi','verbose',0); % set the interal solver to be sedumi
% ops = sdpsettings('solver','sedumi','verbose',0,'sedumi.eps',1e-6);
% ops = sdpsettings('verbose',0,'solver','fmincon'); % nonlinear program
% ops = sdpsettings('solver','sedumi','verbose',0,'sedumi.eps',1e-6,'dualize',1);
 ops = sdpsettings('solver','mosek','cachesolvers',1,'verbose',0,'debug',1,'beeponproblem',1,'saveduals',1); % set the interal solver to be mosek
%  ops = sdpsettings('solver','mosek','cachesolvers',1,'verbose',0,'debug',1,'showprogress',1,'beeponproblem',1,'saveduals',1); % set the interal solver to be mosek
% Call the system parameter function
Systemparameter;
%% If you going to use CVX, then add the following path
% addpath /home/vkien/cvx/cvx/structures;
% addpath /home/vkien/cvx/cvx/lib;
% addpath /home/vkien/cvx/cvx/functions;
% addpath /home/vkien/cvx/cvx/functions/vec;
% addpath /home/vkien/cvx/cvx/commands;
% addpath /home/vkien/cvx/cvx/builtins; 
% addpath /home/vkien/cvx/cvx;
% cvx_setup;

% After delare the matrixes and variables, lets get this party started
% Initiate the random number generators with a random seed
rng('shuffle'); %  the random number generator based on the current time so that rand, randi, and randn produce a different sequence of numbers after each time you call rng.

random_seed = rng; %  seeds the random number generator using the nonnegative integer seed so that rand, randi, and randn produce a predictable sequence of numbers.

Arrival = poissrnd(Mean,K,Iters)'; 
Arrival_Ad = ones(Iters,S);

% Call the MassiveMIMO Only function: where only MBS serves all users
rng(random_seed);
US_Var = ones(1,K+S);
[UserRate_B, Utility_B, Power_Allocate_B ] = MassiveMIMO_Only();

% Call the Scheduling Only function: where MBS  and closed access SCs serve all users
US_Var = ones(1,K+S);
rng(random_seed);
[UserRate_US, Utility_avg_US, Utility_US, Power_Allocate_US ] = SchedulingOnly();

% Call the User Association Only function: where MBS and open access SCs serve all users
rng(random_seed);
% Generate the network localization
%   [MUEnode,  MUEnode_Ad, HeNBnode, HUEnode] = Position( M, S);
  
% Do User Association, Power Allocation, Queue Update

% Save the data and calculate the results based on the trace

% Plot the data

% There is additional function about the Stieltjes transform of nonnegative
% finite measure, Massive MIMO Gain, Beamforming Design
