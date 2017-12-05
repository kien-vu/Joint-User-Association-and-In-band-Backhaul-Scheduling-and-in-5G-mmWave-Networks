function [ UserRate_B, Utility_B, Utility_avg_B, Power_Allocate_B ] = MassiveMIMO_Only()
% Joint Load Balancing and Interference Mitigation in mmWave based Massive MIMO Networks %
% Authors: Trung Kien Vu, 
% Date Update on 2016 March 15
% MATLAB version: 8.4 (R2014b) License host:
%   Username: vkien
%   Host ID: 842b2bc09a54 (eth3)
% CVX Version 2.1, Build 1103 (9714d49)
% Yalmip 
% -----System Parameters---------
% Load system parameter
global ntx % In order to determine the pathloss model

global M  % Number of MUEs
global S  % Number of SCa
global K
global N  % Number of MBS antennas
global F  % Number of antennas is used to mitigate the cross-tier interference
global epsilon % Channel Quality Index
global BW
global Theta_tilde;
global if_fun
global sinr_S2u
global US_Var
global scalingFactor
scalingFactor = 1e6;
global random
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
global Mean % the mean of traffic arrival
global P_b0;
global varphi;
global Tx
global US_Var_A
global N_User
global ops;
global Axis; % In order change the network size 
global Bound1;
Systemparameter;
log_index = 0.69314718056;

Utility_avg_B = zeros(NumTX,N_User);
UserRate_B = zeros(K,Iters,NumTX,N_User);  % average data rare vector of all users
UserS_B = ones(K,Iters,NumTX,N_User); % for users scheduling vector
Power_Allocate_B = zeros(K,Iters,NumTX,N_User); 
%% For Figures
 Rate_avg_B = zeros(Iters,NumTX,N_User);
 Rsue_avg_B = zeros(Iters,NumTX,N_User);
 User_avg_B = zeros(Iters,NumTX,N_User); % MUE 
 Q_avg_B = zeros(Iters,NumTX,N_User);
 Arival_avg_B = zeros(Iters,NumTX,N_User);
 Y_avg_B = zeros(Iters,NumTX,N_User);
 Utility_B =  zeros(Iters,NumTX,N_User);

fprintf('\n Massive MIMO Algorithm is ON\n');


% [MUEnode,  MUEnode_Ad, HeNBnode, HUEnode, THETA, Theta_S, ChanG] = Position( M, S, K, N, F);
for ntx = 1:NumTX
    % Here we increase number of antennas or ratio
    if ntx >= 2
%         M = M + 4;
%         S = S + 4;
% %         K = M + S;
%         if ntx >= 3
%             N = N + 64; % start from N
%         else
%              N = N + 32; % start from N
%         end
%         if ntx >= 7
%              N = N + 100 - 64; % start from N
%         end
    end
%         [MUEnode,  MUEnode_Ad, HeNBnode, HUEnode] = Position_N( M, S);
% Reset the transmit power 
   P_b0 = 26.3794; % 35.9218
for nu = 1:N_User
% Here we increase the distance between MBS and SCs
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
    % Call random seed here to have same setting for each nu. 20160320
%     rng(random_seed);
    disp(nu);


% Variable
r_Scheduling_B = zeros(Iters,K);  % average data rare vector of all users
p_B = zeros(Iters,K);  % average transit power vector to find the auxiliary variable
gamm_1_B = zeros(Iters,K); % auxililary variable Solution 1
GAMA = 4*BW; % set to twice of the mean of traffic arrive
User_Scheduling_B = zeros(Iters,K); % 
% Actual Queues
Q_B = zeros(Iters,K); % Network queue 
% Virtual Queues
Y_B = zeros(Iters,K);  % Set Constraint
YY = zeros(Iters,K);  % Set Constraint for temporary


% Uility function of user capacity 
Rsue = zeros(Iters,S); % The utility function of sue
Omega_B = ones(Iters,K);
A_B = zeros(Iters,K);
CINR_B = zeros(Iters,K);

% This function to call the local spatial correlation matries
% Update on 2016 April
% [THETA, Theta_S, ChanG] = LocalSpatialChannel( S, K, N, F);
% [MUEnode,  MUEnode_Ad, HeNBnode, HUEnode] = Position_N( M, S);


for iter = 1:Iters  
    % Generate the position for each time slot
% [MUEnode,  MUEnode_Ad, HeNBnode, HUEnode] = Position( M, S);
    [MUEnode,  MUEnode_Ad, HeNBnode, HUEnode] = Position_N( M, S);
    [THETA, Theta_S, ChanG] = LocalSpatialChannel( S, K, N, F);


% THETA =  ones(N,N,K); 
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
    Cons1_B = ones(1,K);
    US_Var_A = ones(1,K+S);
    [ ~, THETA_tilde_Only_B, ~, ~, ~, SINR_S2u, ~, ~] = ...
        ChannelModel(MUEnode,  MUEnode_Ad, HeNBnode, HUEnode, THETA, Theta_S, ChanG);
    sinr_S2u = SINR_S2u;
    for k = 1:K
        CINR_B(iter,k) = (1- epsilon.^2);
    end
% We asign the same transmit power to users at the beginning:
    for k = 1:K
        p_B(iter,k) = 10^((P_b0-30)/10)/K;
        UserRate_B(k,iter,ntx,nu) = BW * log( 1 + CINR_B(iter,k) * p_B(iter,k));
        if iter ==1
            gamm_1_B(iter,k) = UserRate_B(k,iter,ntx,nu);
        end
    end
% We model the traffic arrival distribution as an exponential
% distribution with mean mu = 100; exprnd(Mean,K,1); 
% R = poissrnd(lambda) generates random numbers from the Poisson distribution with mean parameter lambda.
% Arrival(1,:) = poissrnd(Mean,K,1)'; 
    for i=1:K
        Q_B(1,i) = max(0 - UserRate_B(i,1,ntx,nu),0) + Arrival(1,i);  % network queue 
        Y_B(1,i) = 0;%+ gamm_1(1,i) - UserRate_B(1,i); % Update queue corresponse to set constraint
    end
    for k=1:K
            r_Scheduling_B(iter,k) = BW * log( 1 +  (CINR_B(iter,k)) * p_B(iter,k));
            User_Scheduling_B(iter,k) = - (Q_B(iter,k) + Y_B(iter,k)) * r_Scheduling_B(iter,k);
    end
        if N >= K
            UserS_B(:,iter,ntx,nu) = ones(1,K)';
            US_Var = ones(1,K+S);
        else
            [US] = UserSchedulingFunction(User_Scheduling_B(iter,:));
            UserS_B(:,iter,ntx,nu) = US';
            US_Var = [US ones(1,S)];
        end        
%         disp('The Massive MIMO Olny user scheduling is FOUND');

%% From now we have a set of user scheduling, we generate the channel model here

%% Calculate the Small user throughput here, since Now we know which SC is
% scheduled
%%  Find the auxiliary variables gamm(iter,i)
        for i = 1:K
            YY(iter, i) =  Y_B(iter, i);
        end
% Solution 1: Using CVX
% Solution 2 More simple things
        for k=1:K
            if Y_B(iter, k) == 0
                gamm_1_B(iter,k) = GAMA;
            else
                gamm_1_B(iter,k) = max(0, min(varphi/(YY(iter, k)), GAMA));
            end
        end      
% disp('The auxiliary variables are FOUND');
%%  Find the data rate vector  R(iter,i) : Interference Mitigation and Power Allocation
    for k=1:K
         A_B(iter,k) = (Q_B(iter,k) * UserS_B(k,iter,ntx,nu) + Y_B(iter,k)) * UserS_B(k,iter,ntx,nu);
    end
%%  Find the Omega vector to calculate the power constraints, this is the Massive MIMO gains
    rng(random1); % call the same channel model as before
    [ THETA_tilde_Wired, THETA_tilde_Only_B, THETA_tilde_Wireless, V, HHH, SINR_S2u, IF_Fun, SINR_S2M] = ...
        ChannelModel(MUEnode,  MUEnode_Ad, HeNBnode, HUEnode, THETA, Theta_S, ChanG);    x0 = ones(1,K);
    Theta_tilde = THETA_tilde_Only_B;
    options = optimset('Display','off');
    [x2,fval,exitflag,output] = fsolve(@OmegaCal,x0,options); % Call solver
    x_B = real(x2);
    for k=1:K
        if x_B(k) >= 1/N
         Omega_B(iter,k) = UserS_B(k,iter,ntx,nu)/(x_B(k));
        else
         Omega_B(iter,k) =  UserS_B(k,iter,ntx,nu) * N;
        end
    end
    
% Scaling down index to [0, 1]
    A_B(iter,:) = real(A_B(iter,:)/(1.5 * max(A_B(iter,:))));
   % Scaling up index to [1, infinity]
   % A_B(iter,:) = A_B(iter,:)/(mean(A_B(iter,:)));
    m = 15; % ensuring that error accurancy is less than 10^-5
% Define variable
    t_p_MIMO       = sdpvar(K,1);
    p_min_MIMO     = sdpvar(K,1);
    kappa_MIMO     = sdpvar(m+4,K,'full');
% Define constraints and Objectives
    constraints = []; % contain all the constraints
    constraints = [constraints, t_p_MIMO >= 0];
    constraints = [constraints, Omega_B(iter,:) * p_min_MIMO - N * 10^((P_b0-30)/10) <= 0];
           for k =1:K
%                 if UserS(k,iter,ntx,nu) == 0
%                     constraints = [constraints, p_min(k) == 0];
%                     constraints = [constraints, t_p(k) == 0];
%                 else
                    constraints = [constraints, p_min_MIMO(k) >= 0];
                    constraints = [constraints, p_min_MIMO(k) <= 10^((P_b0-30)/10)];
                    % Equivalent constrants
%                     constraints = [constraints, CINR(iter,k)* p_min(k) >= exp(t_p(k)) - 1];
                    constraints = [constraints,  kappa_MIMO(:,k) >= 0];
                    constraints = [constraints, 1 + CINR_B(iter,k)* p_min_MIMO(k) >= kappa_MIMO(1,k)];
                    constraints = [constraints, cone([2 + t_p_MIMO(k)/(2^(m-1)); 1-kappa_MIMO(2,k)], kappa_MIMO(2,k) +1)];
                    constraints = [constraints, cone([5/3 + t_p_MIMO(k)/(2^(m)); 1-kappa_MIMO(3,k)], kappa_MIMO(3,k) +1)];
                    constraints = [constraints, cone([2*kappa_MIMO(2,k); 1-kappa_MIMO(4,k)],kappa_MIMO(4,k) +1)];
                    constraints = [constraints, 19/72 + kappa_MIMO(3,k) + 1/24*kappa_MIMO(4,k) <= kappa_MIMO(5,k)];
                    for  mVar = 5:m+3
                        constraints = [constraints, cone([2*kappa_MIMO(mVar,k);1-kappa_MIMO(mVar+1,k)], kappa_MIMO(mVar+1,k) +1)];
                    end
                    constraints = [constraints, cone([2*kappa_MIMO(m+4,k); 1-kappa_MIMO(1,k)], 1+kappa_MIMO(1,k))];
%                 end
           end
     % Define the objective function
     
     obj = - A_B(iter,:) * t_p_MIMO;
    %Solve the problem
    
     sol = optimize(constraints, obj, ops); % solve the problem optimize replaced sdpsolve
         t_p_MIMO = value(t_p_MIMO);
        kappa_MIMO = value(kappa_MIMO);
        p_min_MIMO = value(p_min_MIMO);
     % Check the results
     if sol.problem == 0
         Power_Allocate_B(:,iter,ntx,nu) = value(p_min_MIMO);
     else
         display('PA - MM, something went wrong!');
         sol.info
         yalmiperror(sol.problem)
     end

    clear t_p_MIMO
    clear kappa_MIMO
    clear p_min_MIMO
    clear obj
    clear constraints
    clear diagnostics
     % Solution 2: Using the KKT conditions   
%%  Update the selection of transmit power for next time slot to do user scheduling
    for k=1:K
             UserRate_B(k,iter,ntx,nu) = BW * UserS_B(k,iter,ntx,nu)*  log( 1 + CINR_B(iter,k) * Power_Allocate_B(k,iter,ntx,nu));
             p_B(iter+1,k) = 10^((P_b0-30)/10)/K;%learning * p(iter,k) + (1-learning) * Power_Allocate_B(iter,k); % Using learning method
    end
%%  Update Queues poissrnd
%     Arrival(iter,:) = poissrnd(Mean,K,1)'; % Poission Distribution
%     Arrival(iter,:) = exprnd(Mean,K,1); % Exponential Distribution
%%  Update Queues for the next step

    for k=1:K 
            Q_B(iter+1,k) =  max(Q_B(iter,k) - UserRate_B(k,iter,ntx,nu), 0) + Arrival(iter,k);  % network queue 
            Y_B(iter+1,k) =  max(Y_B(iter,i) + gamm_1_B(iter,k) - UserRate_B(k,iter,ntx,nu), 0); % Update queue corresponse to set constraint
    end
    
end  % iter
 for i=1:Iters
     for k = 1:K
     Rate_avg_B(i,ntx,nu) = Rate_avg_B(i,ntx,nu) + UserRate_B(k,i,ntx,nu);
    
     Q_avg_B(i,ntx,nu) =  Q_avg_B(i,ntx,nu) + Q_B(i,k);
     Arival_avg_B(i,ntx,nu) =  Arival_avg_B(i,ntx,nu) + Arrival(i,k);
     Y_avg_B(i,ntx,nu) = Y_avg_B(i,ntx,nu) + Y_B(i,k);
     end
     for s =1:S
          Rsue_avg_B(i) = Rsue_avg_B(i) + Rsue(i,s);
     end
     for k = 1:K
         if UserRate_B(k,i,ntx,nu) ~= 0
            Utility_B(i,ntx,nu) = Utility_B(i,ntx,nu) +  UserRate_B(k,i,ntx,nu);

         end
     end   
 end
  for i=1:Iters
     for k=1:K
            User_avg_B(i,ntx,nu) = User_avg_B(i,ntx,nu) + UserRate_B(k,i,ntx,nu);
     end
 end
 Rate_avg_B(:,ntx,nu) = Rate_avg_B(:,ntx,nu)/K;
 User_avg_B(:,ntx,nu) = User_avg_B(:,ntx,nu)/K;
 Rsue_avg_B(:,ntx,nu) = Rsue_avg_B(:,ntx,nu)/S;
 Q_avg_B(:,ntx,nu) = Q_avg_B(:,ntx,nu)/K;
 Arival_avg_B(:,ntx,nu) = Arival_avg_B(:,ntx,nu)/K;
 Y_avg_B(:,ntx,nu) = Y_avg_B(:,ntx,nu)/K;
  % Save the immediate results here
for i = 1:Iters
    Utility_avg_B(ntx,nu) = Utility_avg_B(ntx,nu) + Utility_B(i,ntx,nu);
end
    Utility_avg_B(ntx,nu) = Utility_avg_B(ntx,nu)/Iters;
    
%     file0 = strcat('MM_ImmediateResult','_','N',num2str(N),'_','K',num2str(K),'_','atNumTX',num2str(ntx),'_','atN_User_',num2str(nu),'_','P_b0',num2str(P_b0),'_','Tx',num2str(Tx),'BW',num2str(BW),'_','Runs_',num2str(Iters),'On',num2str(date),'.mat');
%     cd Results
%     save(file0,'Utility_avg_B','Utility_B','UserRate_B');
%     cd .. 
end % for increase the  number of users 
end % for increase the number of antennas
 % Save data
  Avg_Q_B = zeros(NumTX,N_User);
 for ntx =1:NumTX
    for nu =1:N_User
        for i=1:Iters
            Avg_Q_B(ntx,nu) = Avg_Q_B(ntx,nu) + Q_avg_B(i,ntx,nu);
        end
            Avg_Q_B(ntx,nu) = Avg_Q_B(ntx,nu)/Iters;
    end
 end 
 
  for ntx =1:NumTX
    for nu =1:N_User
        for i=1:Iters
            Utility_avg_B(ntx,nu) = Utility_avg_B(ntx,nu) + Utility_B(i,ntx,nu);
        end
            Utility_avg_B(ntx,nu) = Utility_avg_B(ntx,nu)/Iters;
    end
  end  
 file_B = strcat('R3Q4MassiveMIMO_Only','N',num2str(N),'_','K',num2str(K),'Ntx',num2str(NumTX),'_','N_User',num2str(N_User),'P_b0',num2str(P_b0),'_','Tx',num2str(Tx),'_','BW',num2str(BW),'_','Runs',num2str(Iters),'_','Random',num2str(rand),'On',num2str(date),'.mat');
    cd Results
    save(file_B);
    cd ..
end

