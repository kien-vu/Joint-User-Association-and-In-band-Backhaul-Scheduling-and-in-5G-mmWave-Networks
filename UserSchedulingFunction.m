function [ US ] = UserSchedulingFunction(User_Scheduling)
%   This function is used to solved the maximum weighted problem. When just
%   select number of users based on the largest weighted factor
    var_K = length(User_Scheduling);
    global N
    US = zeros(1,var_K); % for users scheduling vector% Reset the user scheduling vector at each iter
    for n=1:N
    [Matching, ~] = HungarianMinM(User_Scheduling);
    for k=1:var_K
        if Matching(k) == 1
            US(k) = 1;
            User_Scheduling(k) = 0;
            break;
        end
    end   
    end
end

