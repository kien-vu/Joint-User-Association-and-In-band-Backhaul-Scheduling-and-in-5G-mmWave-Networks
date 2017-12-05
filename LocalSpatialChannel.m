function [THETA, Theta_S, ChanG] = LocalSpatialChannel( S, K, N, F)
%  This function is used to generate the eNodeB and users position.
%  The Covariance matrix also done here
%  Node Position Generation Square Area
% global Axis; % in order change the network size 
% The maximum transmit power should change accordingly to the area. Based
% on the SNR of cell-edge user
% global P_b0;
% global Tx;  % Total transit power at the SC
%%
    % 20 MHz -101 dBm
    % 100 MHz -94.1 dBm
    % 1000 MHz -84 dBm

% We are going to create the covariance matrix*: Antenna arrays
% This is Massive MIMO called: Local Spatial Channel, due to the local
% scateering, the spatial correlation matrices of some users in same are
% mostly identical
THETA =  ones(N,N,K); 
 for k=1:K
    Z = ( randn(N,N) + 1i * randn(N,N))/sqrt(2);
    [R, ~] = qr(Z);
    UN = R(:,1:F); % Semi unitary matrix U' * U = I
    V = ones(1,F);
    B = diag(V); % diagonal matrix D with positive diagonal entries.
    THETA(:,:,k) = UN * B * UN';
 end
 
 Theta_S = ones(N,N,S); 
 for s=1:S
    Z = ( randn(N,N) + 1i * randn(N,N))/sqrt(2);
    [R, ~] = qr(Z);
    UN = R(:,1:F); % Semi unitary matrix U' * U = I
    V = ones(1,F);
    B = diag(V); % diagonal matrix D with positive diagonal entries.
    Theta_S(:,:,s) = UN * B * UN';
 end
  ChanG = ones(1,S);
 for s=1:S
     ChanG(s)  =   sqrt(1/2)* ( randn + 1i*randn);
 end

end