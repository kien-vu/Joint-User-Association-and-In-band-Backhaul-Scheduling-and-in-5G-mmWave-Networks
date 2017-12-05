function [ Opt_US, Opt_US_SUE, Opt_US_MUE, Opt_OP, ConvertP] = UserAssociation(var_UA_MUE, var_SINR_MUE, var_UA_SC, var_SINR_SC, IF_Fun)
% Joint Load Balancing and Interference Mitigation in mmWave based Massive MIMO Networks %
% Authors: Trung Kien Vu, 
% Date Update on 2016 March 15
% MATLAB version: 8.4 (R2014b) License host:
%   Username: vkien
%   Host ID: 842b2bc09a54 (eth3)
% CVX Version 2.1, Build 1103 (9714d49)
% Yalmip 

% Solving user association for 5G ultra-dense small cell networks with wireless backhaul
% The problem is Mixed-Integer Non Convex Program: By using integer relaxation, convex
% approximation and binary search algorithm based on SCA method
% 
% % Varible
global N; % Number of antennas
global K; % Number of MUE and SC users
global F; % Number of data streams
global M; % MUE
global S; % SUE
global N_s; % % Number of MBS antennas of MBS
global N_au % number of active users, SC can serve up to N_au users
global ObjFun; % temporary 
global ops;
global Convergence; % Number of iterations to convert
Convergence = 1;
global ConvP 
global Epsilon; % Interference Threshold
ConvP = 0.001;
L = K + S;
ObjFun = 0;
lambda = ones(L,L); % coefficent between \beta and MUE index
lambda_SU = ones(1,S); % coefficent between \beta and SUE index
lambda_SM = ones(M,S); % coefficent between \beta and MUE index

Cons1 = ones(1,L);
Index = zeros(L,1);
us_min = zeros(L,1); % Temporary variable
% This will get the maximum values of MUEs rates

% Scaling down to [0, 1] and select the best station
ScaleIndex = max(max(max(var_UA_MUE)),max(max(var_UA_SC)));
MacroIndex = var_UA_MUE/(1.5 * ScaleIndex);
SmallIndex = var_UA_SC/( 1.5 * ScaleIndex);
NofIteration = 0; % To count number of iterations

while(1)
     NofIteration = NofIteration + 1;
%      disp(NofIteration)
% Get the FD interference
    for id = 1:L
       if id <= K
            for k = 1:S
                Index(id) =  Index(id) + 1/2 * lambda(id,k) * IF_Fun(id, M + k);
            end
       else
            for k = 1:K
                Index(id) =  Index(id) + 1/(2 * lambda(k,id-K)) * IF_Fun(id - S, k);
            end
       end
    end
% Yalmip solution
% Define variables
    m = 15; % Ensuring that error accurancy is less than 10^-5
    t           =   sdpvar(L,1); % Equivalent variables 
    us_min      =   sdpvar(L,1); % Scheduling variables for MBSs % Note \beta(1:S) = us_min(K+1:L)
    us_min_Mue  =   sdpvar(M,S); % User association variables for MUEs
    t_MUE       =   sdpvar(M,S); % User association variables for MUEs
    us_min_Sue  =   sdpvar(S,1); % Scheduling variables for SUEs
    kappa       =   sdpvar(m+4,L,'full'); % Approximation variables   
    kappaS      =   sdpvar(m+4,S,'full'); % Approximation variables
    kappaM      =   sdpvar(m+4,M,S,'full'); % Approximation variables

   
% Define constraints and Objectives % That fucks me alot
    constraints = []; % contain all the constraints
    obj = - MacroIndex(1,:) * t(1:M) - SmallIndex(1:S)' * t(M+1:K) - ...
          SmallIndex(S+1:2*S)' * ( sqrt(us_min(K+1:L)).* lambda_SU'/2 + 1/2 * sqrt(t(K+1:L)) ./ lambda_SU' ) - ...
          SmallIndex(S+1:2*S)' * ( sqrt(us_min(K+1:L)).* sum(lambda_SM)'/2 + 1/2 * sqrt(sum(t_MUE))' ./ sum(lambda_SM)' ); 
      
% Binary variable constrains
    constraints = [constraints, us_min <= 1];
    constraints = [constraints, us_min >= 0]; % 0
    constraints = [constraints, us_min_Sue <= 1];
    constraints = [constraints, us_min_Sue >= 0];
    constraints = [constraints, us_min_Mue <= 1];
    constraints = [constraints, us_min_Mue >= 0];
    constraints = [constraints, t >= 0];
    constraints = [constraints, t_MUE >= 0];
    
% User scheduling and association constrains for each MUE
    for m = 1:M % Modified each MUE at least will associate with one BS if S + N >= M
       constraints = [constraints, us_min(m) + sum(us_min_Mue(m ,1:S)) <= 1];  % Force the transmission for MUEs
%        constraints = [constraints, us_min(m) + sum(us_min_Mue(m ,1:S)) >= 1];  % Force the transmission for MUEs
       constraints = [constraints, us_min(m) + sum(us_min_Mue(m ,1:S)) >= 0]; 
    end
    
% User scheduling constrains for each SC
    for s = 1:S
       constraints = [constraints, sum( us_min_Mue(:,s)) +  us_min_Sue(s) <= N_au]; % maximum Number of transmissions N_au = 2
       constraints = [constraints, sum( us_min_Mue(:,s)) +  us_min_Sue(s) >= 1]; % Number of transmissions ensuring at least 1 UE to server
       
       if sum( value(us_min_Mue(:,s))) +  value(us_min_Sue(s)) >= 1
         constraints = [constraints, us_min(s + K) == 1]; % Force the transmissions
%          constraints = [constraints, us_min(s + K) -  us_min_Sue(s) == 0]; % Force the transmissions
       end
       constraints = [constraints, us_min(s + K) -  sum( us_min_Mue(:,s)) == 0]; % Force the transmissions
    end   
% Antenna constraints
    constraints = [constraints, sum(Cons1(1:K) * us_min(1:K)) + sum(sum( us_min_Mue)) +  sum(us_min_Sue) - N/F <= 0];
    constraints = [constraints, sum(Cons1(1:K) * us_min(1:K)) + sum(sum( us_min_Mue)) +  sum(us_min_Sue) >= 0];

% User Scheduling Policy - Interference constraints
    constraints = [constraints, cone((us_min .* sqrt(Index)), Epsilon) ];
    
% Equavalent variables
    % For MUEs with MBS index
    m = 15; % Ensuring that error accurancy is less than 10^-5
    for l = 1:M
        constraints = [constraints,  kappa(:,l) >= 0];
        if MacroIndex(l) == 0
            constraints = [constraints, us_min(l) == 0; t(l) == 0];
        else
             constraints = [constraints, 1 + var_SINR_MUE(1,l) * us_min(l)>= kappa(1,l)];
             constraints = [constraints, cone([2 + t(l)/(2^(m-1)); 1-kappa(2,l)], kappa(2,l) +1)];
             constraints = [constraints, cone([5/3 + t(l)/(2^(m)); 1-kappa(3,l)], kappa(3,l) +1)];
             constraints = [constraints, cone([2 * kappa(2,l); 1-kappa(4,l)],kappa(4,l) +1)];
             constraints = [constraints, 19/72 + kappa(3,l) + 1/24*kappa(4,l) <= kappa(5,l)];
                for  mVar = 5:m+3
                    constraints = [constraints, cone([2 * kappa(mVar,l);1 - kappa(mVar+1,l)], kappa(mVar+1,l) +1)];
                end
              constraints = [constraints, cone([2 * kappa(m+4,l); 1 - kappa(1,l)], 1 + kappa(1,l))];

        end
    end
     % For SCs
      for l = M+1:K
        constraints = [constraints,  kappa(:,l) >= 0];
        if SmallIndex(l-M) == 0
            constraints = [constraints, us_min(l) == 0; t(l) == 0];
        else
            constraints = [constraints, 1 + var_SINR_SC(l-M) * us_min(l) >= kappa(1,l)];
            constraints = [constraints, cone([2 + t(l)/(2^(m-1)); 1 - kappa(2,l)], kappa(2,l) +1)];
            constraints = [constraints, cone([5/3 + t(l)/(2^(m)); 1 - kappa(3,l)], kappa(3,l) +1)];
            constraints = [constraints, cone([2 * kappa(2,l); 1 - kappa(4,l)], kappa(4,l) +1)];
            constraints = [constraints, 19/72 + kappa(3,l) + 1/24 * kappa(4,l) <= kappa(5,l)];
            for  mVar = 5:m+3
                constraints = [constraints, cone([2 * kappa(mVar,l);1 - kappa(mVar+1,l)], kappa(mVar+1,l) +1)];
            end
            constraints = [constraints, cone([2 * kappa(m+4,l); 1 - kappa(1,l)], 1 + kappa(1,l))];

        end
       end
    % For SUEs
    for l = 1:S
       constraints = [constraints,  kappaS(:,l) >= 0];
        if SmallIndex(S+l) == 0 % The Queue D
            constraints = [constraints, us_min_Sue(l) == 0; t(l+K) == 0];
        else
             constraints = [constraints, 1 + var_SINR_SC(l+S) * us_min_Sue(l) >= kappaS(1,l)];
             constraints = [constraints, cone([2 + t(l+K)/(2^(m-1)); 1 - kappaS(2,l)], kappaS(2,l) +1)];
             constraints = [constraints, cone([5/3 + t(l+K)/(2^(m)); 1 - kappaS(3,l)], kappaS(3,l) +1)];
             constraints = [constraints, cone([2 * kappaS(2,l); 1 - kappaS(4,l)], kappaS(4,l) +1)];
             constraints = [constraints, 19/72 + kappaS(3,l) + 1/24 * kappaS(4,l) <= kappaS(5,l)];
                for  mVar = 5:m+3
                    constraints = [constraints, cone([2 * kappaS(mVar,l);1 - kappaS(mVar+1,l)], kappaS(mVar+1,l) +1)];
                end
              constraints = [constraints, cone([2 * kappaS(m+4,l); 1 - kappaS(1,l)], 1 + kappaS(1,l))];

        end
    end
     
    % For MUEs with SBSs
    for s = 1:S
        for l = 1:M
            constraints = [constraints,  kappaM(:,l,s) >= 0];
            if SmallIndex(S+s) == 0
                constraints = [constraints, us_min_Mue(l,s) == 0; t_MUE(l,s) == 0];
            else
                 constraints = [constraints, 1 + var_SINR_MUE(s,l) * us_min_Mue(l+s) >= kappaM(1,l,s)];
                 constraints = [constraints, cone([2 + t_MUE(l,s)/(2^(m-1)); 1 - kappaM(2,l,s)], kappaM(2,l,s) +1)];
                 constraints = [constraints, cone([5/3 + t_MUE(l,s)/(2^(m)); 1 - kappaM(3,l,s)], kappaM(3,l,s) +1)];
                 constraints = [constraints, cone([2 * kappaM(2,l,s); 1 - kappaM(4,l,s)], kappaM(4,l,s) +1)];
                 constraints = [constraints, 19/72 + kappaM(3,l,s) + 1/24 * kappaM(4,l,s) <= kappaM(5,l,s)];
                    for  mVar = 5:m+3
                        constraints = [constraints, cone([2 * kappaM(mVar,l,s);1 - kappaM(mVar+1,l,s)], kappaM(mVar+1,l,s) +1)];
                    end
                  constraints = [constraints, cone([2 * kappaM(m+4,l,s); 1 - kappaM(1,l,s)], 1 + kappaM(1,l,s))];

            end
         end
    end
    
% Solve the problem
     sol = optimize(constraints, obj, ops); % solve the problem with optimize command replaced sdpsolve
     us_min = value(us_min);
     t = value(t);
     t_MUE = value(t_MUE);
     us_min_Mue = value(us_min_Mue);
     us_min_Sue = value(us_min_Sue);
     
% Check the results
% Cutting point
cp = 0.001;
if sol.problem == 0
% Update the lambda
   for s = K+1:L
       for k = 1:K
%             if us_min(s) > cp && us_min(k) > cp
%                 lambda(s-K,k) = us_min(k)/us_min(s);
%             elseif us_min(s) <= 0.1 && us_min(k) >= 0.1
%                 lambda(s-K,k) = us_min(k)/2;
%             elseif us_min(k) <= 0.1 && us_min(s) >= 0.1
% %                 lambda(s-K,k) = 2/us_min(s);
%                 lambda(s-K,k) = us_min(s)/2;
%             end
            if us_min(k) <= cp || us_min(s) <= cp
                lambda(s-K,k) = 1;
            else
                lambda(s-K,k) = us_min(k)/us_min(s);
            end
       end
   end
       % Update the lamda carefully. 20160903
   for s = 1:S
       % Update lamda SU
%        if t(M+s) > cp && us_min(M+s) > cp
%            lambda_SU(s) =   us_min(M+s)/t(M+s);
%        elseif  t(M+s) <= 0.1 && us_min(M+s) >= 0.1
%            lambda_SU(s) =   us_min(M+s)/2;
%        elseif t(M+s) >= 0.1 && us_min(M+s) <= 0.1
%            lambda_SU(s) =   t(M+s)/2;
%        end
       if t(M+s) <= cp || us_min(M+s) <= cp
           lambda_SU(s) =  1;
       else
           lambda_SU(s) =   us_min(M+s)/t(M+s);
       end
       % Update lamda MUE
       for m = 1:M
           if  t_MUE(m,s)  <= cp || us_min(M+s) <= cp
               lambda_SM(m,s) =  1;
           else
               lambda_SM(m,s) =   us_min(M+s)/t_MUE(m,s);
           end
       end
    end
% Calculate the objective function
        temp =  MacroIndex(1,1:M) * t(1:M) + SmallIndex(1:S)' * t(M+1:K) + ...
        SmallIndex(S+1:2*S)' * ( sqrt(us_min(K+1:L)).* lambda_SU'/2 + 1/2 * sqrt(t(K+1:L)) ./ lambda_SU' ) + ...
        SmallIndex(S+1:2*S)' * ( sqrt(us_min(K+1:L)).* sum(lambda_SM)'/2 + 1/2 * sqrt(sum(t_MUE))' ./ sum(lambda_SM)' ); 
% Check the convergence
        if( (temp - ObjFun) < ConvP)
%             disp('DONE: The UA Al converges at iteration: ');
%             disp(Convergence)
            break;
        else
%           Increase number of iteration for finding the convergence points
%           disp('Increase iteration to:');
%           disp(Convergence + 1);
            Convergence = Convergence + 1;
            ObjFun = temp;

        end
%          disp('The problem is not accurately solved (or is infeasible) with the current values of lambda');
%          disp('Try to run the algorithm again with different choice of lambda');
%          increase the value of lambda and run again
%          break; % Terminate execution of for or while loop
    else
         display('Hmm, something went wrong!');
         sol.info;
         display(sol.problem);
         yalmiperror(sol.problem);
         us_min = ones(L,1);
         break;
end % end of if condition solved
      
      
end % end of while loop

% Here we deal with integer variables after having the convergence
% point. We numerically observe that the SCA-based Algorithm 1 converges 
% quickly after few iterations and yields a solution of many
% scheduling and operation variables close or equal to binary.
% Hence, we apply a binary search algorithm in order to obtain
% a low-complexity search algorithm to convert the continuous
% relaxation solution to the integer solution


% Opt_US = round(us_min');

% First select which BS will serve MUE
for m =1:M
    if max(us_min_Mue(m,:)) > us_min(m)
        [a,b] = max(us_min_Mue(m,:));
        us_min(m) = 0;
        
        if sum(round(us_min_Mue(1:M,b))) <= 1  % check whether there is another MUE already assigned to SCs b
            us_min_Mue(m,b) = 1;
        else
            us_min_Mue(m,b) = 0;
        end
    else
        us_min(m) = 1;
        for s = 1:S
            us_min_Mue(m,s) = 0;
        end
    end
end
% Second one is to determind the SUEs to be served or not and determine
% beta
us_min_Sue_0 = round(us_min_Sue);
for s = 1:S
    if us_min_Sue_0(s) == 1
        us_min(s+K) = 1;
%     else
%         us_min(s+K) = 0; % maybe remove in case  SCs only serves another MUE expect its users
    end
end
% After select for SCs, then determine the scheduling vector, based on the low search algorithm - Reference below
% H. Li, L. Song, and M. Debbah, “Energy efficiency of large-scale multiple antenna systems
% with transmit antenna selection,” IEEE Transactions on Communications, vol. 62, no. 2, pp. 638–647, 2014.
us_min_0 = us_min';
Certained_Index = zeros(1,L); % return the user index with best response
UnCertained_Index = zeros(1,L); % return the user index with UnCertained response
Removed_Index = zeros(1,L); % return the user index with bad response

% Find the certained users 
i = 1;
jj = 1;
l = 1;
    for  k = 1:L
        if (1.0 - us_min_0(k)) <= 0.25
            Certained_Index(1,i) = k;
            us_min_0(k) = 1;
            i = i+1;
        elseif us_min_0(k) > 0.01 && us_min_0(k) < 0.75
            UnCertained_Index(1,jj) = k;
            jj = jj + 1;
        else
            Removed_Index(1,l) = k;
            us_min_0(k) = 0;
            l = l+1;
        end
    end
% Certained_Index
% UnCertained_Index
% Removed_Index
    J = jj-1; % Number of Uncertained User    
% Rearrange the Ojective function of Uncertained Users in decending order 
    Maximum = zeros(1,J);
    Index = zeros(1,J);
    for j=1:J
        if UnCertained_Index(1,j) <= M
            Maximum(j) = -var_UA_MUE(1,UnCertained_Index(1,j)) * log2( 1+ us_min_0(UnCertained_Index(1,j))*var_SINR_MUE(1,UnCertained_Index(1,j)));
        else
            Maximum(j) = -var_UA_SC(UnCertained_Index(1,j)-M) * log2( 1+ us_min_0(UnCertained_Index(1,j))*var_SINR_SC(UnCertained_Index(1,j)-M));
        end
    end
    for n=1:J
        [Matching, Cost] = HungarianMinM(Maximum); % Need to have - objective function for maximizing the weighted sum function
        for j=1:J
            if Matching(j) == 1
                Index(1,n) = UnCertained_Index(1,j);
                Maximum(1,j) = 0;
                break;
            end
        end   
    end

% Index
for j = 1:J
    % Convert the variable to integer one
    if Index(1,j) <= K % 
        us_min_0(Index(1,j)) = 1;
    % Check the constrains codition to pick user
        Sum1 = 0;
    for ii = 1:i-1
        Sum1 =  Sum1 + us_min_0(Certained_Index(1,ii));
    end
    if (Sum1 + us_min_0(Index(1,j))) <= N
        State = 'true';
        for s = K+1:L
            if (s-K)~= Index(1,j)
                if us_min_0(Index(1,j)) * us_min_0(s) * IF_Fun(s-S,Index(1,j)) - Epsilon > 0
                    State = 'false'; % the user can not satisgy the interference constraint
                end
            end
        end
        if strcmp(State, 'true') == 1
            Certained_Index(1,i) = Index(1,j);
            i = i + 1;
%           disp('Insert uncertained user');
%             disp(Index(1,j));
        else
%             disp('Do not Insert uncertained user');
            Removed_Index(1,l) = Index(1,j);
            us_min_0(Index(1,j)) = 0;
            Index(1,j) = 0; % Remove the uncertained user
            l = l+1;
        end
    else
%           disp('Do not have enough antennas: Break Here');
          break;
    end
    else   % for OP variable
    us_min_0(Index(1,j)) = 1;
% Check the constrains condition to pick user
    Sum1 = 0;
    for ii = 1:i-1
        Sum1 =  Sum1 + us_min_0(Certained_Index(1,ii));
    end
    if (Sum1 + us_min_0(Index(1,j))) <=N
        State = 'true';
        for k=1:K
            if us_min_0(Index(1,j)) * us_min_0(k) * IF_Fun(Index(1,j)-K, k) - Epsilon > 0
                State = 'false'; % the user can not satisgy the interference constraint
            end
        end
        if strcmp(State, 'true') == 1
            Certained_Index(1,i) = Index(1,j);
            i = i + 1;
%           disp('Insert uncertained user');
%           disp(Index(1,j));
        else
%           disp('Do not Insert uncertained user');
            Removed_Index(1,l) = Index(1,j);
            us_min_0(Index(1,j)) = 0;
            Index(1,j) = 0; % Remove the uncertained user
            l = l+1;
        end
    else
%   disp('Do not have enough antennas: Break Here');
    break;
    end

    end % end for index < K

end
% I = i-1; % number of optimal users
% disp('User Scheduling Vector');
% Save the data
% disp('Show the vatiables');
Opt_US_SUE = (us_min_Sue_0)';
Opt_US_MUE = (us_min_Mue)';
Opt_OP = us_min_0(K+1:L);
Opt_US = us_min_0;
ConvertP = Convergence;
end

