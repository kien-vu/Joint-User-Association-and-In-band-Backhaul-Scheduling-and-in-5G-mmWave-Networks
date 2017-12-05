% Joint Load Balancing and Interference Mitigation in mmWave based Massive MIMO Networks %
% Authors: Trung Kien Vu, 
% Date Update on 2016 March 15
% MATLAB version: 8.4 (R2014b) License host:
%   Username: vkien
%   Host ID: 842b2bc09a54 (eth3)
% CVX Version 2.1, Build 1103 (9714d49)
% Yalmip 

% servers
% clc;
% clear;
% close all;
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
% install_sedumi;
% addpath /home/vkien/yalmip/solvers/sdpt3;
% install_sdpt3;
% 
%  addpath /home/vkien/mosek;
%  addpath /home/vkien/mosek/7;
%  addpath /home/vkien/mosek/7/toolbox/r2013a;

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
% rng('shuffle'); %  the random number generator based on the current time so that rand, randi, and randn produce a different sequence of numbers after each time you call rng.

random_seed = rng; %  seeds the random number generator using the nonnegative integer seed so that rand, randi, and randn produce a predictable sequence of numbers.

Arrival = poissrnd(Mean,K,Iters)'; 
Arrival_Ad = ones(Iters,S);
% Initiate the random number generators with a random seed

fprintf('\n Reviewer 3 Question 3\n');

rng(random_seed);
US_Var = ones(1,K+S);
[UserRate_B, Utility_B, Utility_avg_B, Power_Allocate_B ] = MassiveMIMO_Only();

Utility_avg_B
US_Var = ones(1,K+S);
rng(random_seed);
[UserRate_US, Utility_avg_US, Utility_US, Power_Allocate_US ] = SchedulingOnly();


rng(random_seed);
disp('random number');
disp(rng);
Systemparameter;
disp('User Association Algorithm is ON');
for ntx=1:NumTX
    % Here we increase number of antennas or ratio
%     if ntx >= 2
% %         M = M + 4;
% %         S = S + 4;
% %         K = M + S;
%         if ntx >= 3
%             N = N + 64; % start from N
%         else
%              N = N + 32; % start from N
%         end
%         if ntx >= 7
%              N = N + 100 - 64; % start from N
%         end
%     end
%     [MUEnode,  MUEnode_Ad, HeNBnode, HUEnode] = Position_N( M, S);
% Reset the transmit power 
    P_b0 = 41; 
for nu = 1:N_User
    
%        P_b0 = 15 + PathLoss(Axis) - 84; % in dB
    if nu > 1
%          Var = Var + 5; % start from N
%            BW = BW + BW *(ntx-1);
%         Bound2 = Bound2 + 150;
%       M = M + 2;
        P_b0 = P_b0 - 3; % in dB
% Note if you decrease the area, you need to scale the transmit power too. 
%         Axis = Axis - 100;
%         P_b0 = 15 + PathLoss(Axis/2) - 84; % 84 in dB at 28GHz
%         Bound1 = Bound1 - 5;
    end
    disp('Run time at nu');
    disp(nu);
    
    % Call random seed here to have same setting for each nu. 20160320
    rng(random_seed);


P_k = 10^((P_b0-30)/10)/K; % Total transit power at the macro BS
P_s = 10^((Tx-30)/10); % Total transit power at the macro BS

SINR = zeros(1,Iters);
SINR_DE = zeros(1,Iters);
% Systemparameter;



p = zeros(Iters,K);  % average transit power vector to find the auxiliary variable
p_A = zeros(Iters,K);  % average transit power vector to find the auxiliary variable
gamm_1 = zeros(Iters,K); % auxililary variable Solution 1
gamm_A = zeros(Iters,K); % auxililary variable Solution 2
GAMA = 4 * BW; % set to twice of the mean of traffic arrive

User_Association = zeros(Iters,K+S); % user for User scheduling, outside index of log function
SINR_Scheduling = zeros(Iters,K+S);  % user for User scheduling, inside index of log function
User_Association_HD = zeros(Iters,K); % user for User scheduling, outside index of log function
SINR_Scheduling_HD = zeros(Iters,K);  % user for User scheduling, inside index of log function
% Actual Queues
Q = zeros(Iters,K); % Network queue 
% Virtual Queues
Y = zeros(Iters,K);  % Virtual queue
YY = zeros(Iters,K);  % Set Constraint for temporary
D = zeros(Iters,S);  % Wireless backhaul constraint: SC queue
  
% Uility function of user capacity 
R = ones(Iters,K); % The utility function of all users: MUEs and SCs
g = zeros(Iters,1); % The function of power constraint
Omega = ones(Iters,K);
x = ones(Iters,K);
x_1 = ones(Iters,M); % just find for M users for HetNets Wired case
B = ones(1,K);
A = zeros(Iters,K);
CINR = zeros(Iters,K);
Op_Tx_Sc = zeros(Iters,S);

for iter = 1:Iters  
    % Generate the position for each time slot
    [MUEnode,  MUEnode_Ad, HeNBnode, HUEnode] = Position_N( M, S);
    [THETA, Theta_S, ChanG] = LocalSpatialChannel( S, K, N, F);
     random1 = rng;
     US_Var = ones(1,K+S);
     US_Var_A = ones(1,K+S);
    [THETA_tilde_Wired, THETA_tilde_Only_B, THETA_tilde_Wireless, V, HHH, SINR_S2u, IF_Fun, SINR_S2M] = ...
        ChannelModel(MUEnode,  MUEnode_Ad, HeNBnode, HUEnode, THETA, Theta_S, ChanG);
    if_fun = IF_Fun;
    sinr_S2u = SINR_S2u;
%     SINR_S2M
    % Precalculate the data rate for small cell user: taking account for FD
    % interference
    for s = 1:S
        Rsue(s,iter,ntx,nu)  =  BW  * log(1 + N_au * SINR_S2u(s)/(1 + Epsilon)); % poissrnd(Mean);  %
    end 
    % Channel Quality
    CINR(iter,1:M) = (1- epsilon.^2);
%     CINR(iter,M+1:K) = 10^(Tx_Gain/10);
    CINR(iter,M+1:K) = 1;

% We asign the same transmit power to users at the beginning:
    for k=1:K
        p(iter,k) = 10^((P_b0-30)/10)/K;%
        UserRate(k,iter,ntx,nu) = BW * log2( 1 + CINR(iter,k) * p(iter,k) );
        if k > M
           UserRate(k,iter,ntx,nu) = BW * log2( 1 + CINR(iter,k) * p(iter,k) * 10^(Tx_Gain/10) );
        end
        if iter ==1
            gamm_1(k,iter,ntx,nu) = UserRate(k,iter,ntx,nu);
        end
    end
    % Queue Initilization
    for i=1:K
        Q(1,i) = max(1 - UserRate(i,1,ntx,nu),1) + Arrival(1,i);  % network queue 
        Y(1,i) = 1; % Update virtual queue corresponse to set constraint
    end
    for s=1:S % Update the SC queue at the first time slot
        D(1,s) = max(1 +  UserRate(s+M,1,ntx,nu) - Rsue(s,1,ntx,nu),1) + Arrival(1,s+M);% Update queue corresponse to wireless backhaul queue constraint
    end 
    
    % Finding variable for User Scheduling and Operation mode
    L = K + S;
    for k = 1:L
        if k <= M
            SINR_Scheduling_MUE(1,k,iter,ntx,nu)   =   CINR(iter,k) * p(iter,k) /(1 + S*0.005);
            User_Association_MUE(1,k,iter,ntx,nu)  =   ( Q(iter,k)  + Y(iter,k) );
            for b = 2:(1+S)
                 SINR_Scheduling_MUE(b,k,iter,ntx,nu)   =  SINR_S2M(k,b-1);
                 User_Association_MUE(b,k,iter,ntx,nu)  = ( Q(iter,k)  + Y(iter,k) );
            end
        elseif M < k && k <= K
            SINR_Scheduling_SC(k-M,iter,ntx,nu)   =  CINR(iter,k) * p(iter,k);
            User_Association_SC(k-M,iter,ntx,nu)  = ( Q(iter,k)  + Y(iter,k))  /(1 + S*0.005) ;
        else
            SINR_Scheduling_SC(k-M,iter,ntx,nu)   =   SINR_S2u(k-K) /(1 + (S-1)*0.005); % Note to change the formulation SINR_S2u = p |h|^2
            User_Association_SC(k-M,iter,ntx,nu)  =  D(iter, k-K);
        end
    end
    
% Call the User Association for Load Balancing + OP function : Update on 2016 March 7
    [ Opt_US, Opt_US_SUE, Opt_US_MUE, Opt_OP, ConvertP] = UserAssociation( User_Association_MUE(:,1:M,iter,ntx,nu), ...
        SINR_Scheduling_MUE(:,1:M,iter,ntx,nu), User_Association_SC(:,iter,ntx,nu), SINR_Scheduling_SC(:,iter,ntx,nu), if_fun);
    UserS(:,iter,ntx,nu) = Opt_US';
    UA_MUE(:,:,iter,ntx,nu) = Opt_US_MUE;
% Update the convergence point of SCA method
    ConvergencePoints(iter,ntx,nu) = ConvertP;
    NuofMUEs_associatedwithSCs(iter,ntx,nu) = sum(sum(Opt_US_MUE));



%% From now we have a set of user scheduling, we generate the channel model here
%% Calculate the Small user throughput here, since Now we know which SC is
% scheduled
for s=1:S
    Rsue(s,iter,ntx,nu) =  UserS(K + s, iter,ntx,nu) * BW  * log(1 + (SINR_S2u(s) * 10^(Tx_Gain/10) )/(1 + Epsilon)); %poissrnd(Mean);  %
end 
   

%%  Find the auxiliary variables gamm(iter,i)
for i = 1:K
    if i > M
        YY(iter, i) =  Y(iter, i) + D(iter, i- M);
    else
        YY(iter, i) =  Y(iter, i);
    end
end

for k = 1:K
        if Y(iter, k) == 0
            gamm_1(iter,k) = GAMA;
        else
            if k > M
                temp1 = 0; % Total data rate arrive for MUEs
                for m = 1:M
                    temp1 = temp1 + BW * Opt_US_MUE(k-M,m)  *  log( 1 + SINR_Scheduling_MUE(k-M,m,iter,ntx,nu) * P_s/2);
                end 
                gamm_1(iter,k) = max(Rsue(k-M,iter,ntx,nu) + temp1, min( varphi/(YY(iter, k) ), GAMA));
            else
                gamm_1(iter,k) = max(0, min( varphi/(YY(iter, k) ), GAMA));
            end
        end
end       
%     disp('The auxiliary variables are FOUND');
    
%%  Find the data rate vector  R(iter,i) : Interference Mitigation and Power Allocation
    for k = 1:K
         A(iter,k) = ( Q(iter,k) * UserS(k,iter,ntx,nu) + Y(iter,k) * UserS(k,iter,ntx) );
    end
    
%%  Find the Omega vector to calculate the power constraints 
    rng(random1); % call the same channel model as before, we have to call the channel model function again, ...
    % since now we selected a new subset of users to design the precoder and assign the transmit power
    US_Var = UserS(:,iter,ntx,nu);
if sum(US_Var(1:K)) ~= 0 % Check whether having users to serve   
    [ THETA_tilde_Wired, THETA_tilde_Only_B, THETA_tilde_Wireless, V, HHH, SINR_S2u, IF_Fun, SINR_S2M] = ...
    ChannelModel(MUEnode,  MUEnode_Ad, HeNBnode, HUEnode, THETA, Theta_S, ChanG);
    x0 = ones(1,K);
    options = optimset('Display','off');
    Theta_tilde = THETA_tilde_Wireless;
    [x_C,fval,exitflag,output] = fsolve(@OmegaCal,x0,options); % Call solver
    x(iter,:) = real(x_C);
    for k=1:K
        if x(iter,k) >= 1/N
         Omega(iter,k) = UserS(k,iter,ntx,nu)/x(iter,k);
        else
         Omega(iter,k) =  UserS(k,iter,ntx,nu) * N;
        end
    end
% Find Power Allocation Vector  
% Scaling index down [0, 1] to ensure the convex
    A(iter,:) = real(A(iter,:)/(1.5 * max(A(iter,:))));
% Using CVX solution
% Yalmip solution using SOCP
    m = 15; % ensuring that error accurancy is less than 10^-5
% Define variable
    t_p       = sdpvar(K,1);
    p_min     = sdpvar(K,1);
    kappa     = sdpvar(m+4,K,'full');
% Define constraints and Objectives
    constraints = []; % contain all the constraints
    constraints = [constraints, t_p >= 0];
    constraints = [constraints, Omega(iter,:) * p_min - N * 10^((P_b0-30)/10) <= 0];
           for k = 1:K
                if UserS(k,iter,ntx,nu) == 0
                    constraints = [constraints, p_min(k) == 0];
                    constraints = [constraints, t_p(k) == 0];
                else
                    constraints = [constraints, p_min(k) >= 0];
                    constraints = [constraints, p_min(k) <= 10^((P_b0-30)/10)];
                % Equivalent constrants
                    constraints = [constraints,  kappa(:,k) >= 0];
                    constraints = [constraints, 1 + CINR(iter,k)* p_min(k) >= kappa(1,k)];
                    constraints = [constraints, cone([2 + t_p(k)/(2^(m-1)); 1-kappa(2,k)], kappa(2,k) +1)];
                    constraints = [constraints, cone([5/3 + t_p(k)/(2^(m)); 1-kappa(3,k)], kappa(3,k) +1)];
                    constraints = [constraints, cone([2*kappa(2,k); 1-kappa(4,k)],kappa(4,k) +1)];
                    constraints = [constraints, 19/72 + kappa(3,k) + 1/24*kappa(4,k) <= kappa(5,k)];
                    for  mVar = 6:m+4
                        constraints = [constraints, cone([2*kappa(mVar-1,k);1-kappa(mVar,k)], kappa(mVar,k) +1)];
                    end
                    constraints = [constraints, cone([2*kappa(m+4,k); 1-kappa(1,k)], 1+kappa(1,k))];
                end
           end
     obj = - A(iter,:) * t_p;
% Solve the problem
     sol = optimize(constraints, obj, ops); % solve the problem optimize replaced sdpsolve
     p_min = value(p_min);
% Check the results
     sol.info;
     if sol.problem == 0
         Power_Allocate(:,iter,ntx,nu) = (p_min);
     else
         display('PA, something went wrong!');
         sol.info;
         yalmiperror(sol.problem)
     end

else
        Power_Allocate(:,iter,ntx,nu) = 0;
end
   clear kappa
   clear p_min
   clear obj
   clear constraints
   clear diagnostics
         
%%  Update the selection of transmit power for next time slot to do user association 20163015
%   Opt_US, Opt_US_SUE, Opt_US_MUE, Opt_OP, ConvertP
    for k = 1:K
        if k <= M
            UserRate(k,iter,ntx,nu) = BW * UserS(k,iter,ntx,nu) *  log( 1 + CINR(iter,k) * Power_Allocate(k,iter,ntx,nu));
            temp = 0;
            for s = 1:S
                temp = max(temp,  BW * Opt_US_MUE(s,k)  *  log( 1 + SINR_Scheduling_MUE(s,k,iter,ntx,nu) * P_s/2) );
            end
            UserRate(k,iter,ntx,nu) = max (temp, UserRate(k,iter,ntx,nu));
        else
            UserRate(k,iter,ntx,nu) = BW * UserS(k,iter,ntx,nu)*  log( 1 + CINR(iter,k) * Power_Allocate(k,iter,ntx,nu) * 10^(Tx_Gain/10) );
        end
            p(iter+1,k) = 10^((P_b0-30)/10)/K; % learning * p(iter,k) + (1-learning) * Power_Allocate(iter,k); % Using learning method
    end

%%  Update Queues for the next step
    for k = 1:K   
        Q(iter+1,k) = max(Q(iter,k) - UserRate(k,iter,ntx,nu), 1) + Arrival(iter,k);  % network queue 
        Y(iter+1,k) = max(Y(iter,k) + gamm_1(iter,k) - UserRate(k,iter,ntx,nu), 1); % Update queue corresponse to set constraint
    end
    % New update here: More UEs per SC
    for s = 1:S
        temp = 0; % Total data rate arrive for MUEs
        for m = 1:M
            temp = temp + BW * Opt_US_MUE(s,m)  *  log( 1 + SINR_Scheduling_MUE(s,m,iter,ntx,nu) * P_s/2);
        end
        D(iter + 1,s) =  max(D(iter,s) + UserRate(M+s,iter,ntx,nu) - ( Rsue(s,iter,ntx,nu) + temp), 1); % Update queue corresponse to wireless backhaul constraint
    end 
end  % end of iter
 
 for i = 1:Iters
     for k = 1:K
         Rate_avg(i,ntx,nu) = Rate_avg(i,ntx,nu) + UserRate(k,i,ntx,nu);
         Q_avg(i,ntx,nu) =  Q_avg(i,ntx,nu) + Q(i,k);
         Arival_avg(i,ntx,nu) =  Arival_avg(i,ntx,nu) + Arrival(i,k);
         Y_avg(i,ntx,nu) = Y_avg(i,ntx,nu) + Y(i,k);
     end
     for s = 1:S
        D_avg(i,ntx,nu) =  D_avg(i,ntx,nu) + D(i,s);
        Rsue_avg(i,ntx,nu) = Rsue_avg(i,ntx,nu) + Rsue(s,i,ntx,nu);
        R_SC(i, ntx,nu) = R_SC(i, ntx,nu) + UserRate(s+M,i,ntx,nu);

     end
     for k = 1:K
         if  UserRate(k,i,ntx,nu) > 0
            Utility(i,ntx,nu) = Utility(i,ntx,nu) +  UserRate(k,i,ntx,nu);
         end
     end  
     
 end
for i = 1:Iters
     for k = 1:K
        User_avg(i,ntx,nu) = User_avg(i,ntx,nu) +  UserRate(k,i,ntx,nu);
        if k <= M
            User_avg(i,ntx,nu) = User_avg(i,ntx,nu) + UserRate(k,i,ntx,nu);
        else
            User_avg(i,ntx,nu) = User_avg(i,ntx,nu) + min(Rsue(k-M,i,ntx,nu), UserRate(k,i,ntx,nu));
        end
     end
     for m = 1:M
         R_MUE(i, ntx,nu) = R_MUE(i, ntx,nu) + UserRate(m,i,ntx,nu);
     end
end
    

% Rate Average
 Rate_avg(:,ntx,nu) = Rate_avg(:,ntx,nu);
 User_avg(:,ntx,nu) = User_avg(:,ntx,nu);
 R_MUE(:,ntx,nu) = R_MUE(:, ntx,nu);
 Rsue_avg(:,ntx,nu) = Rsue_avg(:,ntx,nu);
 % Queue Length Avearge
 Q_avg(:,ntx,nu) = Q_avg(:,ntx,nu)/K; % calculate the aggregate queue: average for each user.
 Arival_avg(:,ntx,nu) = Arival_avg(:,ntx,nu);
 Y_avg(:,ntx,nu) = Y_avg(:,ntx,nu)/K; % calculate the aggregate queue: average for each user.
 D_avg(:,ntx,nu) = D_avg(:,ntx,nu)/S; % calculate the aggregate queue: average for each user.

 % Save the immediate results here
for i = 1:Iters
    Utility_avg_UA(ntx,nu) = Utility_avg_UA(ntx,nu) + Utility(i,ntx,nu);
end
    Utility_avg_UA(ntx,nu) = Utility_avg_UA(ntx,nu)/Iters;
%     file0 = strcat('UA_ImmediateResult','_','N',num2str(N),'_','K',num2str(K),'_','atNumTX',num2str(ntx),'_','atN_User_',num2str(nu),'_','P_b0',num2str(P_b0),'_','Tx',num2str(Tx),'BW',num2str(BW),'_','Runs_',num2str(Iters),'_','Random',num2str(rand),'On',num2str(date),'.mat');
%     cd Results
%     save(file0,'Utility_avg_UA','Utility','UserRate');
%     cd ..
end % for increase the number of users
end % for increase the number of antennas



 for ntx = 1:NumTX
    for nu = 1:N_User
        for i = 1:Iters
            Utility_avg_UA(ntx,nu) = Utility_avg_UA(ntx,nu) + Utility(i,ntx,nu);
        end
            Utility_avg_UA(ntx,nu) = Utility_avg_UA(ntx,nu)/Iters;
    end
 end   
 Avg_Q = zeros(NumTX,N_User);
 for ntx = 1:NumTX
    for nu = 1:N_User
        for i = 1:Iters
            Avg_Q(ntx,nu) = Avg_Q(ntx,nu) + Q_avg(i,ntx,nu);
        end
            Avg_Q(ntx,nu) = Avg_Q(ntx,nu)/Iters;
    end
 end 

   Avg_D_UA = zeros(NumTX,N_User);
 for ntx=1:NumTX
    for nu =1:N_User
        for i=1:Iters
        Avg_D_UA(ntx,nu) = Avg_D_UA(ntx,nu) + D_avg(i,ntx,nu);
        end
        Avg_D_UA(ntx,nu) = Avg_D_UA(ntx,nu)/Iters;
    end
 end 


% Save the data
file1 = strcat('R3Q4','_','N',num2str(N),'_','K',num2str(K),'_','NumTX_',num2str(NumTX), ...
    '_','NUser_',num2str(N_User),'_','P_b0',num2str(P_b0),'_', ...
    'Tx_',num2str(Tx),'Axis_',num2str(Axis),'_','Runs_',num2str(Iters),'_','Random_',num2str(rand),'On',num2str(date),'.mat');
cd Results
save(file1);
cd ..

fprintf('\n Finished \n');

