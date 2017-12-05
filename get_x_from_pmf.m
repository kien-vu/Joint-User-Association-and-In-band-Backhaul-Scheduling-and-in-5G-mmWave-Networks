function n = get_x_from_pmf(x, pmf)
% this function is going to pick one action from the set of actions
% according to the probability distribution (t-1), eq(16)
% x - numbers
% pmf - probability mass function
% n - random number

r = rand(1);% test random
cdf = 0;
if ( round(sum(pmf)) == 1) && (length(x) == length(pmf) )
    for i = 1:length(pmf)
        cdf = cdf + pmf(i);
        if (r-cdf<=0) %random value is between pmf(i-1) and pmf(i)
            break;
        end
    end
    n = x(i);
else
    n = null;
end