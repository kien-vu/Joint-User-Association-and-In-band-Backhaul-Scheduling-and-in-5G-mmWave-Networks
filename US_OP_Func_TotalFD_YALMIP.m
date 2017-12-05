function [ Opt_US ] = US_OP_Func_TotalFD_YALMIP(var_US, var_OM, IF_Fun)
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
ops = sdpsettings('verbose',0,'solver','fmincon'); % nonlinear program
% define output
Opt_US = zeros(L,1);
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
%            Index
%          cvx_begin quiet
% %        cvx_solver gurobi 
%         cvx_solver sedumi %sedumi sdpt3
%         variable us_min(L)
%         variable t(L)
%             maximize geo_mean(t);
%             %maximize var_US*log ( 1 + var_OM' .* us_min )
%             subject to
% % Binary variable constrains
%                  us_min <= 1;
%                  us_min >= 0;
% % Equivalent constrants
%                  for l=1:L
%                      if var_US(l) == 0
%                          us_min(l) == 0;
%                          t(l) == 0;
%                      else
%                      var_OM(l) * us_min(l) >= power(t(l), (1/var_US(l))) - 1;
%                      end
%                  end
% % Explicit constrant
%                  t >= 0;
% % Antenna constraints
%            sum(Cons1*us_min) - N/F <= 0;
% % User Scheduling Policy % Total FD interference constraint
%             sum(square(us_min) .* Index) - Epsilon <=0;
%         cvx_end 
        
% Yalmip solution
% Define variable
    t       = sdpvar(L,1);
    us_min  = sdpvar(L,1);
% Define constraints and Objectives
    constraints = []; % contain all the constraints
    constraints = [constraints, t >= 0];
    obj = -geomean(t);%- var_US*slog ( 1 + var_OM' .* us_min );% - slog ( 1 + var_OM' .* us_min );  % -geo_mean(t,1); %- var_US*slog ( 1 + var_OM' .* us_min ); 
% Binary variable constrains
    constraints = [constraints, us_min <= 1];
    constraints = [constraints, us_min >= 0];
    constraints = [constraints, t >= 0];
% Antenna constraints
    constraints = [constraints, sum(Cons1*us_min) - N <= 0];
% User Scheduling Policy % Total FD interference constraint % square by
% norm
    constraints = [constraints, sum(norm(us_min) .* Index) - Epsilon <=0 ];
% Equivalent constrants
    for l = 1:L
        if var_US(l) == 0
            constraints = [constraints, us_min(l) == 0; t(l) == 0];
        else
            constraints = [constraints, power(t(l), (1/var_US(l))) - 1 <= var_OM(l) * us_min(l)];
        end
    end
         diagnostics = solvesdp(constraints, obj, ops); % solve the problem   
         us_min = value(us_min);
%          diagnostics.info
        switch diagnostics.problem
            case 0 % successfully solved
          % Update the lambda
           for k=1:K
               for s=1:S
                    if us_min(k) ~= 0 && us_min(K+s) ~= 0
                        lambda(k,s) = us_min(K+s)/us_min(k);
                    elseif us_min(k) == 0 && us_min(K+s) ~= 0
                        lambda(k,s) = us_min(K+s)/2;
                    elseif us_min(K+s) == 0 && us_min(k) ~= 0
                        lambda(k,s) = 2/us_min(k);
                    end
                    if us_min(K+s) == 0 && us_min(k) == 0
                        lambda(k,s) = 1;
                    end
               end
           end
            % Calculate the objective function
            temp = sum(var_US*log ( 1 + var_OM' .* us_min ));
%             ObjFun
%             us_min
            if(abs(temp-ObjFun) < ConvP)
                disp('FD-Algorithm converges at iteration: ');
                disp(Convergence)
                break;
            else
            % Increase number of iteration for finding the convergence
            % points
                Convergence = Convergence + 1;
                ObjFun = temp;
                
            end
            case 1 % infeasible problem

                break;
                
            case 4 % numerical problem in the internal solver
                disp('There is a numerical problem in the internal solver');
                disp('The algorithm teminates before reaching the maximum number of iterations');
                break;
                
            otherwise % unknown problem
                disp('Failed. Unknown problem -- Unbounded objective function');
                return;
        end
    
 %         if(strfind(cvx_status,'Solved'))
%           % Update the lambda
%            for k=1:K
%                for s=1:S
%                     if us_min(k) ~= 0 && us_min(K+s) ~= 0
%                         lambda(k,s) = us_min(K+s)/us_min(k);
%                     elseif us_min(k) == 0 && us_min(K+s) ~= 0
%                         lambda(k,s) = us_min(K+s)/2;
%                     elseif us_min(K+s) == 0 && us_min(k) ~= 0
%                         lambda(k,s) = 2/us_min(k);
%                     end
%                     if us_min(K+s) == 0 && us_min(k) == 0
%                         lambda(k,s) = 1;
%                     end
%                end
%            end
%             % Calculate the objective function
%             temp = sum(var_US*log ( 1 + var_OM' .* us_min ));
% %             ObjFun
% %             us_min
%             if(abs(temp-ObjFun) < ConvP)
% %                 disp('FD-Algorithm converges at iteration: ');
% %                 disp(Convergence - 1)
%                 break;
%             else
%             % Increase number of iteration for finding the convergence
%             % points
%                 Convergence = Convergence + 1;
%                 ObjFun = temp;
%                 
%             end
% %         else   
% %             disp('The problem is not accurately solved (or is infeasible) with the current values of lambda');
% %             disp('Try to run the algorithm again with different choice of lambda');
% %             % increase the value of lambda and run again
% %             break; % Terminate execution of for or while loop
%         
%         end
        
 end % end while

        % Here we deal with integer variables after having the convergence
        % point. By using the Branch and Bound Method
        % Second option is to use the Binary Search Algorithm
us_min_0 = us_min';
Certained_Index = zeros(1,N); % return the user index with best response
UnCertained_Index = zeros(1,N); % return the user index with UnCertained response
Removed_Index = zeros(1,N); % return the user index with bad response
% Find the certained users 
i = 1;
j = 1;
l = 1;
    for  k = 1:length(us_min_0)
        if (1.0 - us_min_0(k)) < 0.5
            Certained_Index(1,i) = k;
            us_min_0(k) = 1;
            i = i+1;
        elseif us_min_0(k) > 0.01 && us_min_0(k) < 0.5
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
    Maximum(j) = us_min_0(UnCertained_Index(1,j))*var_US(UnCertained_Index(1,j))*log2( 1+ var_OM(UnCertained_Index(1,j)));
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
                    Sum1 =  Sum1 + Cons1(Certained_Index(1,ii))*us_min_0(Certained_Index(1,ii));
                    %Sum1 =  Sum1 + us_min_0(Certained_Index(1,ii))*var_US(Certained_Index(1,ii))*log2( 1+ var_OM(Certained_Index(1,ii)));
                end
%                 Sum1
                if (Sum1 + Cons1(Index(1,j))*us_min_0(Index(1,j))) <=N
                    State = 'true';
                    for s=K+1:L
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
for id =1:I
    Opt_US(Certained_Index(id)) = 1;
end



end

