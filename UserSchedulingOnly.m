function [ UserRate_US, Utility_avg_US, Utility_US, Power_Allocate_US ] = UserSchedulingOnly()

% Joint Self-Backhauling and Interference Mitigation in Massive MIMO 5G Networks%
% Authors: Trung Kien Vu, 
% Date Update on 2015 May 06
% Date Update on 2015 Set 11
% Date Update on 2016 March 12, Saturday morning
% MATLAB version: 8.4 (R2014b) License host:
%   Username: vkien
%   Host ID: 842b2bc09a54 (eth3)
% CVX Version 2.1, Build 1103 (9714d49)
% Yalmip too

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
global scalingFactor
scalingFactor = 1e6;
% random = rng;% used for generate the same channel parameter for comaprision
% Initial Value
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
global random_seed
global Var
% Load CVX command
% Initial Value
global Iters; % or time slot
global NumTX;
global P_b0;
global varphi;
global Tx
global US_Var_A
global N_User
global Bound2
global ops;
global if_fun
global Tx_Gain
global N_s 
global N_au 
omega = 1; % weight of proportionally fair utility function
log_index = 0.69314718056;
global Mean
Systemparameter;

%% For Figures
 Rate_avg_US = zeros(Iters,NumTX,N_User);
 Rsue_avg_US = zeros(Iters,NumTX,N_User);
 User_avg_US = zeros(Iters,NumTX,N_User); % MUE and SUE (min(SC,SUE))
 Q_avg_US = zeros(Iters,NumTX,N_User);
 Arival_avg_US = zeros(Iters,NumTX,N_User);
 Y_avg_US = zeros(Iters,NumTX,N_User);
 D_avg_US = zeros(Iters,NumTX,N_User);
 Utility_US =  zeros(Iters,NumTX,N_User);
 R_MUE_US = zeros(Iters,NumTX,N_User);
 R_SC = zeros(Iters,NumTX,N_User);


UserS_US = ones(K+S,Iters,NumTX,N_User); % for users scheduling vector
Power_Allocate_US = zeros(K,Iters,NumTX,N_User); % Using CVX
UserRate_US = zeros(K,Iters,NumTX,N_User); % Using CVX  % average data rare vector of all users
Rsue_US = zeros(S,Iters,NumTX,N_User); % The utility function of sue


[MUEnode,  MUEnode_Ad, HeNBnode, HUEnode, THETA, Theta_S, ChanG] = Position( M, S, K, N, F);
for ntx = 1:NumTX
    % Here we increase number of antennas or ratio
    if ntx >= 2
%         M = M + 4;
%         S = S + 4;
%         K = M + S;
        if ntx >= 5
            N = N + 32; % start from N
        else
             N = N + 16; % start from N
        end
        if ntx > 10
             N = N + 64 - 32; % start from N
        end
    end
for nu = 1:N_User
% Here we increase the distance between MBS and SCs
    if nu > 1
%          Var = Var + 5; % start from N
%          BW = BW + BW *(ntx-1);
%         Bound2 = Bound2 + 150;
        P_b0 = P_b0 - 3; % in dB
        disp('Run time at nu');
        disp(nu);
    end
    % Call random seed here to have same setting for each nu. 20160320
%     rng(random_seed);
P_k = 10^((P_b0-30)/10)/K; % Total transit power at the macro BS
P_s = 10^((Tx-30)/10); % Total transit power at the macro BS

p_US = zeros(Iters,K);  % average transit power vector to find the auxiliary variable
gamm_1_US = zeros(Iters,K); % auxililary variable Solution 1
GAMA_US = 4 * BW; % set to twice of the mean of traffic arrive

User_Scheduling_US = zeros(Iters,K+S); % user for User scheduling, outside index of log function
OP_Scheduling_US = zeros(Iters,K+S);  % user for User scheduling, inside index of log function

% Actual Queues
Q_US = zeros(Iters,K); % Network queue 
% Virtual Queues
Y_US = zeros(Iters,K);  % Virtual queue
YY_US = zeros(Iters,K);  % Set Constraint for temporary
D_US = zeros(Iters,S);  % Wireless backhaul constraint: SC queue

% Uility function of user capacity 
Omega_US = ones(Iters,K);
x_US = ones(Iters,K);

A_US = zeros(Iters,K);
CINR_US = zeros(Iters,K);
% We model the traffic arrival distribution as an exponential
% distribution with mean mu = 100; exprnd(Mean,K,1); 
% R = poissrnd(lambda) generates random numbers from the Poisson distribution with mean parameter lambda.
% Arrival(1,:) = poissrnd(Mean,K,1)'; 
% Arrival_Ad(1,:) = poissrnd(Mean,S,1)'; 
% [MUEnode,  MUEnode_Ad, HeNBnode, HUEnode, THETA, Theta_S, ChanG] = Position( M, S, K, N, F);

for iter = 1:Iters  
    % Generate the position for each time slot
    [MUEnode,  MUEnode_Ad, HeNBnode, HUEnode, THETA, Theta_S, ChanG] = Position( M, S, K, N, F);
%     THETA =  ones(N,N,K); 
%  for k=1:K
%     Z = ( randn(N,N) + 1i * randn(N,N))/sqrt(2);
%     [R, P] = qr(Z);
%     U = R(:,1:F); % Semi unitary matrix U' * U = I
%     V = ones(1,F);
%     B = diag(V); % diagonal matrix D with positive diagonal entries.
%     THETA(:,:,k) = U * B * U';
%  end
     random1 = rng;
     US_Var = ones(1,K+S);
     US_Var_A = ones(1,K+S);
     [ THETA_tilde_Wired, THETA_tilde_Only_B, THETA_tilde_Wireless, V, HHH, SINR_S2u, IF_Fun, SINR_S2M] = ...
        ChannelModel(MUEnode,  MUEnode_Ad, HeNBnode, HUEnode, THETA, Theta_S, ChanG);
    sinr_S2u = SINR_S2u;
    if_fun = IF_Fun;
    % Precalculate the data rate for small cell user: taking account for FD
    % interference
    for s = 1:S % times by N_au since there is only one SUE
        Rsue_US(s,iter,ntx,nu)  =  BW  * log(1 + N_au * SINR_S2u(s) /(1 + Epsilon)); % poissrnd(Mean);  %
    end 
    % Channel Quality
    CINR_US(iter,1:M) = (1- epsilon.^2);
    CINR_US(iter,M+1:K) = 10^(Tx_Gain/10);

% We asign the same transmit power to users at the beginning:
    for k = 1:K
        p_US(iter,k) = 10^((P_b0-30)/10);%
        UserRate_US(k,iter,ntx,nu) = BW * log( 1 + CINR_US(iter,k) * p_US(iter,k) );
        if iter == 1
            gamm_1_US(iter,k) = UserRate_US(k,iter,ntx,nu);
        end
    end
   
    % Queue Initilization
    for i=1:K
        Q_US(1,i) = max(0 - UserRate_US(i,1,ntx,nu),1) + Arrival(1,i);  % network queue 
        Y_US(1,i) = 0; % Update virtual queue corresponse to set constraint
    end
    for s=1:S % Update the SC queue at the first time slot
        D_US(1,s) = max(0 +  UserRate_US(s+M,1,ntx,nu) - Rsue_US(s,1,ntx,nu),1) + Arrival(1,s+M);% Update queue corresponse to wireless backhaul queue constraint
    end 
    % Finding variable for User Scheduling and Operation mode
    L = K + S;
    for k = 1:L
        if k <= M
            OP_Scheduling_US(iter,k) =(CINR_US(iter,k)) * p_US(iter,k)/(1 + Epsilon);
            User_Scheduling_US(iter,k) = ( Q_US(iter,k) + Y_US(iter,k) ) ;
        elseif M < k && k <= K
            OP_Scheduling_US(iter,k) = (CINR_US(iter,k)) * p_US(iter,k)/(1 + Epsilon);
            User_Scheduling_US(iter,k) = ( Q_US(iter,k) + Y_US(iter,k)) ;
        else
            OP_Scheduling_US(iter,k) = SINR_S2u(k-K) /(1 + Epsilon); % Note to change the formulation SINR_S2u = p |h|^2
            User_Scheduling_US(iter,k) = D_US(iter,k-K);
        end
    end
    % Call the US+ OP function: Proposal
%     if (K+S) > N
        [ Opt_US_US ] = US_OP_Func_TotalFD_YALMIP_approx(User_Scheduling_US(iter,:), OP_Scheduling_US(iter,:), if_fun);
%       [ Opt_US1 ] = US_OP_Func_TotalFD_SOCP(User_Scheduling(iter,:), OP_Scheduling(iter,:),if_fun);
         UserS_US(:,iter,ntx,nu) = Opt_US_US';
%         Opt_US'
%     else
%         UserS_US(:,iter,ntx,nu) = ones(K+S,1);
%     end
    

%% From now we have a set of user scheduling, we generate the channel model here
%% Calculate the Small user throughput here, since Now we know which SC is
% scheduled
    for s=1:S
        Rsue_US(s,iter,ntx,nu) =  UserS_US(K +s,iter,ntx,nu) * BW  * log(1 + (SINR_S2u(s) * 10^(Tx_Gain/10) )/(1 + Epsilon)); %poissrnd(Mean);  %
    end 
    
%%  Find the auxiliary variables gamm(iter,i)
        for i = 1:K
            if i > M
                YY_US(iter, i) = Y_US(iter, i) + D_US(iter, i- M);
            else
                YY_US(iter, i) =  Y_US(iter, i);
            end
        end

        for k=1:K
                if Y_US(iter, k) == 0
                    gamm_1_US(iter,k) = GAMA_US;
                else
                    gamm_1_US(iter,k) = max(0, min(varphi/(YY_US(iter, k) * log(2) ), GAMA_US));
                end
        end       
%     disp('The auxiliary variables are FOUND');

%%  Find the data rate vector  R(iter,i) : Interference Mitigation and Power Allocation
    for k=1:K
         A_US(iter,k) = ( Q_US(iter,k)*UserS_US(k,iter,ntx,nu) + Y_US(iter,k)*UserS_US(k,iter,ntx) ) / log_index;
    end
%%  Find the Omega vector to calculate the power constraints 

    rng(random1); % call the same channel model as before
    US_Var = UserS_US(:,iter,ntx,nu);
if sum(US_Var(1:K)) ~=0 % Check whether having users to serve   
    [ THETA_tilde_Wired, THETA_tilde_Only_B, THETA_tilde_Wireless, V, HHH, SINR_S2u, IF_Fun, SINR_S2M] = ...
        ChannelModel(MUEnode,  MUEnode_Ad, HeNBnode, HUEnode, THETA, Theta_S, ChanG);
    x0 = ones(1,K);
    options = optimset('Display','off');
    Theta_tilde = THETA_tilde_Wireless;
    [x_C,fval,exitflag,output] = fsolve(@OmegaCal,x0,options); % Call solver
    x_US(iter,:) = real(x_C);
    for k=1:K
        if x_US(iter,k) >= 1/N
         Omega_US(iter,k) = UserS_US(k,iter,ntx,nu)/x_US(iter,k);
        else
         Omega_US(iter,k) =  UserS_US(k,iter,ntx,nu) * N;
        end
    end
%   Find Power Allocation Vector  
% Scaling index down [0, 1] to ensure the convex
    A_US(iter,:) = real(A_US(iter,:)/(1.5 * max(A_US(iter,:))));
%   Yalmip solution using SOCP
    m = 10; % ensuring that error accurancy is less than 10^-5
%   Define variable
    t_us       = sdpvar(K,1);
    p_min_us     = sdpvar(K,1);
    kappa_us     = sdpvar(m+4,K,'full');
% Define constraints and Objectives
    constraints_us = []; % contain all the constraints
    constraints_us = [constraints_us, t_us >= 0];
    constraints_us = [constraints_us, Omega_US(iter,:) * p_min_us - N * P_b0 <= 0];
   for k =1:K
       if  UserS_US(k,iter,ntx,nu) == 0
            constraints_us = [constraints_us, p_min_us(k) == 0];
            constraints_us = [constraints_us, t_us(k)   == 0];
       else
            constraints_us = [constraints_us, p_min_us(k) >= 0];
            constraints_us = [constraints_us, p_min_us(k) <= P_b0];
            % Equivalent constrants
            constraints_us = [constraints_us,  kappa_us(:,k) >= 0];
            constraints_us = [constraints_us, 1 + CINR_US(iter,k)* p_min_us(k) >= kappa_us(1,k)];
            constraints_us = [constraints_us, cone([2 + t_us(k)/(2^(m-1)); 1-kappa_us(2,k)], kappa_us(2,k) +1)];
            constraints_us = [constraints_us, cone([5/3 + t_us(k)/(2^(m)); 1-kappa_us(3,k)], kappa_us(3,k) +1)];
            constraints_us = [constraints_us, cone([2*kappa_us(2,k); 1-kappa_us(4,k)],kappa_us(4,k) +1)];
            constraints_us = [constraints_us, 19/72 + kappa_us(3,k) + 1/24*kappa_us(4,k) <= kappa_us(5,k)];
            for  mVar = 5:m+3
                constraints_us = [constraints_us, cone([2*kappa_us(mVar,k);1-kappa_us(mVar+1,k)], kappa_us(mVar+1,k) +1)];
            end
            constraints_us = [constraints_us, cone([2*kappa_us(m+4,k); 1-kappa_us(1,k)], 1+kappa_us(1,k))];
       end
   end
     obj = - A_US(iter,:) * t_us;
     %Solve the problem
     sol = optimize(constraints_us, obj, ops); % solve the problem optimize replaced sdpsolve
     % Check the results
%      check(constraints)
     if sol.problem == 0
         Power_Allocate_US(:,iter,ntx,nu) = value(p_min_us);
%          disp('power allocation solved');
     else
         display('US of PA something went wrong!');
         sol.info;
         yalmiperror(sol.problem);
     end

    else
        Power_Allocate_US(:,iter,ntx,nu) = 0;
end
   p_min_us = value(p_min_us);
   t_us = value(t_us);
   kappa_us = value(kappa_us);
   clear t_us
   clear kappa_us
   clear p_min_us
   clear obj
   clear constraints_us
   clear diagnostics 
%%  Update the selection of transmit power for next time slot to do user scheduling
    for k=1:K
            UserRate_US(k,iter,ntx,nu) = BW * UserS_US(k,iter,ntx,nu)*  log( 1 + CINR_US(iter,k) * Power_Allocate_US(k,iter,ntx,nu));
%             UserRate(k,iter,ntx,nu) = BW * UserS(k,iter,ntx,nu)*  log2( 1 + CINR(iter,k) * Power_Allocate(k,iter,ntx,nu)/NN0/(Sum_Interference_FD(k)/NN0 + 1/NN0));
            if k > M
              UserRate_US(k,iter,ntx,nu) = BW * UserS_US(k,iter,ntx,nu)*  log( 1 + (CINR_US(iter,k) * Power_Allocate_US(k,iter,ntx,nu)) );
%                UserRate(k,iter,ntx,nu) = BW * UserS(k,iter,ntx,nu)*  log2( 1 + CINR(iter,k) * Power_Allocate(k,iter,ntx,nu));
            end
            p_US(iter+1,k) = 10^((P_b0-30)/10)/K; % learning * p(iter,k) + (1-learning) * Power_Allocate(iter,k); % Using learning method
    end

%%  Update Queues for the next step
    for i=1:K   
            Q_US(iter+1,i) = max(Q_US(iter,i) - UserRate_US(i,iter,ntx,nu), 1) + Arrival(iter,i);  % network queue 
            Y_US(iter+1,i) = max(Y_US(iter,i) + gamm_1_US(iter,i) - UserRate_US(i,iter,ntx,nu), 1); % Update queue corresponse to set constraint
    end

    for s=1:S % Here is only one SUE, then we have the aggegate queue at SC as below
        D_US(iter+1,s) =  max(D_US(iter,s) + UserRate_US(M+s,iter,ntx,nu) - Rsue_US(s,iter,ntx,nu), 1); % Update queue corresponse to wireless backhaul constraint
    end 
end  % end of iter
 
 for i = 1:Iters
     for k = 1:K
     Rate_avg_US(i,ntx,nu) = Rate_avg_US(i,ntx,nu) + UserRate_US(k,i,ntx,nu);
     Q_avg_US(i,ntx,nu) =  Q_avg_US(i,ntx,nu) + Q_US(i,k);
     Arival_avg_US(i,ntx,nu) =  Arival_avg_US(i,ntx,nu) + Arrival(i,k);
     Y_avg_US(i,ntx,nu) = Y_avg_US(i,ntx,nu) + Y_US(i,k);
     end
     for s = 1:S
          D_avg_US(i,ntx,nu) =  D_avg_US(i,ntx,nu) + D_US(i,s);
          Rsue_avg_US(i,ntx,nu) = Rsue_avg_US(i,ntx,nu) + Rsue_US(s,i,ntx,nu);
          R_SC(i, ntx,nu) = R_SC(i, ntx,nu) + UserRate_US(s+M,i,ntx,nu);

     end
     for k = 1:K
         if  UserRate_US(k,i,ntx,nu) > 0
            Utility_US(i,ntx,nu) = Utility_US(i,ntx,nu) +  UserRate_US(k,i,ntx,nu);
         end
     end  
     
 end
for i = 1:Iters
     for k=1:K
        if k <= M
            User_avg_US(i,ntx,nu) = User_avg_US(i,ntx,nu) + UserRate_US(k,i,ntx,nu);
        else
            User_avg_US(i,ntx,nu) = User_avg_US(i,ntx,nu) + min(Rsue_US(k-M,i,ntx,nu), UserRate_US(k,i,ntx,nu));

        end
     end
     for m=1:M
         R_MUE_US(i, ntx,nu) = R_MUE_US(i, ntx,nu) + UserRate_US(m,i,ntx,nu);
     end
end

 Rate_avg_US(:,ntx,nu) = Rate_avg_US(:,ntx,nu);
 User_avg_US(:,ntx,nu) = User_avg_US(:,ntx,nu);
 R_MUE_US(:, ntx,nu) = R_MUE_US(:, ntx,nu);
 
 Rsue_avg_US(:,ntx,nu) = Rsue_avg_US(:,ntx,nu);
 Q_avg_US(:,ntx,nu) = Q_avg_US(:,ntx,nu)/K; % calculate the aggregate queue: average for each user.
 Arival_avg_US(:,ntx,nu) = Arival_avg_US(:,ntx,nu);
 Y_avg_US(:,ntx,nu) = Y_avg_US(:,ntx,nu)/K; % calculate the aggregate queue: average for each user.
 D_avg_US(:,ntx,nu) = D_avg_US(:,ntx,nu)/S; % calculate the aggregate queue: average for each user.
 
end % for increase the number of users
end % for increase the number of antennas

% Calculate real time average data rate for SUEs = min(time average of DL  MBS-SC, time average of DL  SC-SUE )
 Arg_SUE = zeros(NumTX,N_User);
 for ntx=1:NumTX
    for nu =1:N_User       
        Arg_SUE(ntx,nu) = min(sum( Rsue_avg_US(:,ntx,nu)),sum( R_SC(:,ntx,nu)));  
        Arg_SUE(ntx,nu) = Arg_SUE(ntx,nu)/Iters;
    end
 end   

 Utility_avg_US = zeros(NumTX,N_User);
 for ntx=1:NumTX
    for nu =1:N_User
        for i=1:Iters
        Utility_avg_US(ntx,nu) = Utility_avg_US(ntx,nu) + Utility_US(i,ntx,nu);
        end
        Utility_avg_US(ntx,nu) = Utility_avg_US(ntx,nu)/Iters;
    end
 end   


 Avg_Q_US = zeros(NumTX,N_User);
 for ntx=1:NumTX
    for nu =1:N_User
        for i=1:Iters
        Avg_Q_US(ntx,nu) = Avg_Q_US(ntx,nu) + Q_avg_US(i,ntx,nu);
        end
        Avg_Q_US(ntx,nu) = Avg_Q_US(ntx,nu)/Iters;
    end
 end 
  Avg_D_US = zeros(NumTX,N_User);
 for ntx=1:NumTX
    for nu =1:N_User
        for i=1:Iters
        Avg_D_US(ntx,nu) = Avg_D_US(ntx,nu) + D_avg_US(i,ntx,nu);
        end
        Avg_D_US(ntx,nu) = Avg_D_US(ntx,nu)/Iters;
    end
 end 

file3 = strcat('SchedulingOnlyData','_','N',num2str(N),'_','K',num2str(K),'_','NumTX',num2str(NumTX),'P_b0',num2str(P_b0),'_','Tx',num2str(Tx),'_','N_User_',num2str(N_User),'_','BW',num2str(BW),'_','Runs_',num2str(Iters),'On',num2str(date),'.mat');
save(file3);
end


 