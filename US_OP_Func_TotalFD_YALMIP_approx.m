function [ Opt_US ] = US_OP_Func_TotalFD_YALMIP_approx(var_US, var_OM, IF_Fun)
% User scheduling and operation mode for conference version
%   Detailed explanation goes here
% % Authors: Trung Kien Vu (vkien@ee.oulu.fi)
% Address: Center for Wireless Communucations
% Date Update on 2015 May 06
% MATLAB version: 8.4 (R2014b) License host:
%   Username: vkien
%   Host ID: 842b2bc09a54 (eth3)
% Using YALMIP tool, version 20141218

% Solving User Scheduling and operation mode of SCs for full duplex mode
% The problem is Mixed-Integer Non Convex Program: By using integer relaxation, success convex
% approximation and binary search algorithm
% % Varible
global N; % Number of antennas
global K; % Number of MUE and SC users
global F; % Number of data streams
global M;
global S;
global L;
L = K + S;
global ConvP 
global ops
ConvP = 0.005;
lambda = ones(L,L);
Cons1 = ones(1,L);
Index = zeros(L,1);
global Epsilon; % Interference Threshold
global ObjFun;
ObjFun = 0;
global Convergence
Convergence = 1;
% Scaling index down [0, 1] to ensure the convex
  var_US = var_US/(1.5 *max(var_US));
% ops = sdpsettings('solver','sdpt3','verbose',0); % set the interal solver to be SDPT3
% ops = sdpsettings('solver','sedumi','verbose',0); % set the interal solver to be sedumi
% ops = sdpsettings('solver','mosek','verbose',0); % set the interal solver to be mosek
% ops = sdpsettings('verbose',0,'solver','fmincon'); % nonlinear program
% define output
  while(1)
           for id = 1:L
               if id <= K
                    for k = 1:S
                        Index(id) =  Index(id) + 1/2 * lambda(id,k)*IF_Fun(id,M+k);
                    end
               else
                    for k = 1:K
                        Index(id) =  Index(id) + 1/(2 * lambda(k,id-K))*IF_Fun(id - S,k);
                    end
               end
           end   
% Yalmip solution
% Define variable
    m = 15; % ensuring that error accurancy is less than 10^-5
    t_US       = sdpvar(L,1);
    us_min_US  = sdpvar(L,1);
    kappa   = sdpvar(m+4,L,'full');
% Define constraints and Objectives
    constraints = []; % contain all the constraints
    obj = - var_US * t_US; 
% Binary variable constrains
    constraints = [constraints, us_min_US <= 1];
    constraints = [constraints, us_min_US >= 0];
    constraints = [constraints, t_US >= 0];
% Antenna constraints
    constraints = [constraints, sum(Cons1 * us_min_US) - N <= 0];
% User Scheduling Policy % Total FD interference constraint % square by norm
     constraints = [constraints, sum(norm(us_min_US) .* Index) - Epsilon <=0 ];
% Equivalent constrants
    m = 15;
    for l = 1:L
           constraints = [constraints,  kappa(:,l) >= 0];
        if var_US(l) == 0
            constraints = [constraints, us_min_US(l) == 0; t_US(l) == 0];
        else
             constraints = [constraints, 1 + var_OM(l) * us_min_US(l)>= kappa(1,l)];
             constraints = [constraints, cone([2 + t_US(l)/(2^(m-1)); 1-kappa(2,l)], kappa(2,l) +1)];
             constraints = [constraints, cone([5/3 + t_US(l)/(2^(m)); 1-kappa(3,l)], kappa(3,l) +1)];
             constraints = [constraints, cone([2 * kappa(2,l); 1 - kappa(4,l)], kappa(4,l) +1)];
             constraints = [constraints, 19/72 + kappa(3,l) + 1/24 * kappa(4,l) <= kappa(5,l)];
                for  mVar = 5:m+3
                    constraints = [constraints, cone([2 * kappa(mVar,l); 1-kappa(mVar+1,l)], kappa(mVar+1,l) +1)];
                end
              constraints = [constraints, cone([2 * kappa(m+4,l); 1 - kappa(1,l)], 1 + kappa(1,l))];

        end
    end
     %Solve the problem
     sol = optimize(constraints, obj, ops); % solve the problem with optimize command replaced sdpsolve
     % Check the results
     % Cutting point
       cp = 0.001;
     if sol.problem == 0
            us_min_US = value(us_min_US);
            t_US = value(t_US);
          % Update the lambda
           for k=1:K
               for s=1:S
                    if us_min_US(k) <= cp || us_min_US(K+s) <= cp
                        lambda(k,s) = 1;
                    elseif us_min_US(k) == cp && us_min_US(K+s) ~= cp
                        lambda(k,s) = us_min_US(K+s)/2;
                    elseif us_min_US(K+s) == cp && us_min_US(k) ~= cp
                        lambda(k,s) = 2/us_min_US(k);
                    else
                        lambda(k,s) = us_min_US(K+s)/us_min_US(k);
                    end
               end
           end
            % Calculate the objective function 
            temp = - value(obj);
            if( (temp - ObjFun) < ConvP)
%                 disp('US-Hybird-Algorithm converges at iteration: ');
%                 disp(Convergence)
                break;
            else
            % Increase number of iteration for finding the convergence
            % points
                Convergence = Convergence + 1;
                ObjFun = temp; 
            end
    else
         display('FD for Scheduling, something happens!');
         sol.info
         yalmiperror(sol.problem)
    end  
 end % end while


        % Here we deal with integer variables after having the convergence
        % point. By using the Branch and Bound Method
        % Second option is to use the Binary Search Algorithm
us_min_0 = us_min_US';
 clear us_min_US
 clear t_US
 
Certained_Index = zeros(1,N); % return the user index with best response
UnCertained_Index = zeros(1,N); % return the user index with UnCertained response
Removed_Index = zeros(1,N); % return the user index with bad response
% Find the certained users 
i = 1;
j = 1;
l = 1;
    for  k = 1:length(us_min_0)
        if (1.0 - us_min_0(k)) <= 0.25
            Certained_Index(1,i) = k;
            us_min_0(k) = 1;
            i = i+1;
        elseif us_min_0(k) > 0.01 && us_min_0(k) < 0.75
            UnCertained_Index(1,j) = k;
            j = j+1;
        else
            Removed_Index(1,l) = k;
            us_min_0(k) = 0;
            l = l+1;
        end
    end
%     Certained_Index
%     UnCertained_Index
J = j-1; % Number of Uncertained User    
% Rearrange the Ojective function of Uncertained Users in decending order 
Maximum = zeros(1,J);
Index = zeros(1,J);
%minimize var_US*log ( 1 + var_OM' .* us_min )
for j=1:J
    Maximum(j) = - us_min_0(UnCertained_Index(1,j)) * var_US(UnCertained_Index(1,j)) * log( 1+ var_OM(UnCertained_Index(1,j)));
end
for n=1:J
    [Matching, Cost] = HungarianMinM(Maximum);
    for j=1:J
        if Matching(j) == 1
%             Del=k;
            Index(1,n) = UnCertained_Index(1,j);
            Maximum(1,j) = 0;
            break;
        end
    end   
end
for j=1:J
    % Convert the variable to integer one
                if Index(1,j) <= K % 
                us_min_0(Index(1,j)) = 1;
                % Check the constrains codition to pick user
                Sum1 = 0;
                for ii = 1:i-1
                    Sum1 =  Sum1 + Cons1(Certained_Index(1,ii)) * us_min_0(Certained_Index(1,ii));
                    %Sum1 =  Sum1 + us_min_0(Certained_Index(1,ii))*var_US(Certained_Index(1,ii))*log2( 1+ var_OM(Certained_Index(1,ii)));
                end
%                 Sum1
                if (Sum1 + Cons1(Index(1,j)) * us_min_0(Index(1,j))) <= N
                    State = 'true';
                    for s = K+1:L
                        if us_min_0(Index(1,j))*us_min_0(s) * IF_Fun(s-S,Index(1,j)) - Epsilon > 0
                            State = 'false'; % the user can not satisgy the constraint
                        end
                    end
%                     disp(false);
                    if strcmp(State, 'true') == 1
                        Certained_Index(1,i) = Index(1,j);
                        i = i + 1;
%                         disp('Insert uncertained user');
%                         disp(Index(1,j));
                    else
%                         disp('Do not Insert uncertained user');
                        Removed_Index(1,l) = Index(1,j);
                        us_min_0(Index(1,j)) = 0;
                        Index(1,j) = 0; % Remove the uncertained user
                        l = l+1;
                    end
                else
%                     disp('Do not have enough antennas: Break Here');
%                     break;
                end
                else   % for OP variable
                    
                    
                us_min_0(Index(1,j)) = 1;
                % Check the constrains condition to pick user
                Sum1 = 0;
                for ii = 1:i-1
                    Sum1 =  Sum1 + Cons1(Certained_Index(1,ii))*us_min_0(Certained_Index(1,ii));
                    %Sum1 =  Sum1 + us_min_0(Certained_Index(1,ii))*var_US(Certained_Index(1,ii))*log2( 1+ var_OM(Certained_Index(1,ii)));
                end
                if (Sum1 + Cons1(Index(1,j))*us_min_0(Index(1,j))) <=N
                    State = 'true';
                    for id=1:K
                        if us_min_0(Index(1,j))*us_min_0(id) * IF_Fun(Index(1,j)-K, id) - Epsilon > 0
                            State = 'false'; % the user can not satisgy the constraint
                        end
                    end
                    if strcmp(State, 'true') ~= 1
                        Certained_Index(1,i) = Index(1,j);
                        i = i + 1;
%                         disp('Insert uncertained user');
%                         disp(Index(1,j));
                    else
%                         disp('Do not Insert uncertained user');
                        Removed_Index(1,l) = Index(1,j);
                        us_min_0(Index(1,j)) = 0;
                        Index(1,j) = 0; % Remove the uncertained user
                        l = l+1;
                    end
                else
%                     disp('Do not have enough antennas: Break Here');
%                     break;
                end
                    
                end % end for index < K

end
I = i-1; % number of optimal users
% disp('User Scheduling Vector');
% Certained_Index
Opt_US = zeros(1,L);
for id = 1:I
    Opt_US(Certained_Index(id)) = 1;
end
end

