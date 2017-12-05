function [Fun,fun] = OmegaCal(x)
% Calculate the Stieltjes transform of a nonnegative finite measure: Omega
% For each time cal the function, THETA donot allow to change: HOW (global variable)
Size = length(x);
Fun = zeros(Size,1);
global Theta_tilde;
global N;
THETA = Theta_tilde;
SF_0 = max(max(max(abs(THETA))));
THETA = THETA/SF_0;

alpha = 0.01; % This is the RZF parameter 0.01
fun = zeros(N,N);
for k=1:Size
    fun(:,:) = fun(:,:) +  THETA(:,:,k)/(alpha/SF_0 + x(k)/SF_0);
end
for k =1:Size
    Fun(k) =   x(k) - (1/N*trace( THETA(:,:,k) * ( fun/N + eye(N) )^(-1) ));
end
% We solve the equation here