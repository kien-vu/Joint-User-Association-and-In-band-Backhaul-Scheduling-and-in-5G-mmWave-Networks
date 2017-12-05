function [ THETA_tilde_Wired, THETA_tilde_Only_B, THETA_tilde_Wireless, V, HHH, SINR_S2u, IF_Fun, SINR_S2M] = ChannelModel(MUEnode,  MUEnode_Ad, HeNBnode, HUEnode, THETA, Theta_S, ChanG)
% Channel generation
global Axis; % in order change the network size 
global P_b0;
global Tx  % Total transit power at the SC
X_Axis = Axis; % the size of network area
Y_Axis = Axis;
global rounding
rounding = 12;
global K
global M
global S
global F
global N
global US_Var
global US_Var_A
global epsilon
global Tx_Gain
global N_s % % Number of MBS antennas of MBS
global N_au % number of active users, SC can serve up to N_au users
global Epsilon; % Interference Threshold

P_k = 10^((P_b0-30)/10)/K; % Total transit power at the macro BS
P_s = 10 ^((Tx-30)/10); % Total transit power at the macro BS
MeNB = [X_Axis/2; Y_Axis/2];
pathloss_UT = zeros(1,K); % Including MUEs and SCs
pathloss_UT_dB = zeros(1,K);
pathloss_UT1 = zeros(1,S); 
pathloss_UT1_dB = zeros(1,S);
pathloss_UT_Ad_dB = zeros(1,S);
pathloss_UT_Ad = zeros(1,S); % Including new MUEs
 HM = zeros(N,M);  % Channel matrix between MeNB and MUE
 HS = zeros(N,S);  % Channel matrix between MeNB and SC
 HC = zeros(N,S);  % Channel matrix between MeNB and SUE
UT = zeros(2,K);
for k=1:M
    UT(:,k) = MUEnode(:,k);
    pathloss_UT_dB(k) = PathLoss(norm(MeNB - UT(:,k)));
    pathloss_UT(k) = 10^(-pathloss_UT_dB(k)/10);
end
for k=(M+1):(M+S)
    UT(:,k) = HeNBnode(:,k-M);
    pathloss_UT_dB(k) = PathLoss(norm(MeNB - UT(:,k)));
    pathloss_UT(k) = 10^(-pathloss_UT_dB(k)/10);
    pathloss_UT1_dB(k-M) = PathLoss(norm(MeNB - UT(:,k)));
    pathloss_UT1(k-M)=  10^(-pathloss_UT1_dB(k-M)/10);
end
for s=1:S
     pathloss_UT_Ad_dB(k) = PathLoss(norm(MeNB - MUEnode_Ad(:,s)));

    pathloss_UT_Ad(s) = 10^(-pathloss_UT_Ad_dB(s)/10);
end

%% From HERE to calculate Normalized spatial correlation matrics : Solve it
THETA_tilde_Only = ones(N,N,K); % the channel correlation matrix of MUE and SC users 
THETA_tilde_Only_A = ones(N,N,K); % the channel correlation matrix of MUE and SC users 
THETA_tilde_Only_B = ones(N,N,K); % the channel correlation matrix of MUE and SC users 
 for k=1:K
    THETA_tilde_Only(:,:,k) = US_Var(k) * pathloss_UT(k) *  THETA(:,:,k) / trace( THETA(:,:,k)); %pathloss_UT(k) *
    if k <= M
        THETA_tilde_Only_B(:,:,k) = US_Var(k) * pathloss_UT(k)  *  THETA(:,:,k) / trace( THETA(:,:,k)); %same as old MUE
        THETA_tilde_Only_A(:,:,k) = US_Var_A(k) * pathloss_UT(k)  *  THETA(:,:,k) / trace( THETA(:,:,k)); %same as old MUE
    else
        THETA_tilde_Only_B(:,:,k) = US_Var(k) * pathloss_UT1(k-M) *  THETA(:,:,k) / trace( THETA(:,:,k)); %same as old MUE
        THETA_tilde_Only_A(:,:,k) = US_Var_A(k) * pathloss_UT1(k-M) *  THETA(:,:,k) / trace( THETA(:,:,k)); % like SC without antenna gain
    end
 end
 
  for n=1:N
     for m=1:M
         HM(n,m) = sqrt(N) * (sqrt(1-epsilon^2)*(0/N*randn + 1i*1/N*randn) + epsilon*(0/N*randn + 1i*1/N*randn));
     end
     for s=1:S
         HS(n,s) =  sqrt(N) * ( 0*randn + 1i*1/N*randn);
         HC(n,s) =  sqrt(N) * ( 0*randn + 1i*1/N*randn);
     end
 end
 
 
 for m=1:M
     HM(:,m) = sqrt(THETA_tilde_Only(:,:,m)) * HM(:,m);
 end
 for s=1:S
     HS(:,s) = sqrt( THETA_tilde_Only(:,:,M+s)) * HS(:,s);
 end
 THETA_MUE1 = ones(N,N,S); % the channel correlation matrix of additional MUEs only for the HetNet Wired Scenario
 for k=1:S
    THETA_MUE1(:,:,k) = US_Var_A(k+M) * pathloss_UT_Ad(k) * N * Theta_S(:,:,k) / trace(Theta_S(:,:,k)); %pathloss_UT(k) *
 end
 THETA_MUE_Ad = ones(N,N,M+S); % the channel correlation matrix of MUEs only for the HetNet Wired Scenario
 THETA_MUE_Ad(:,:,1:K) = THETA_tilde_Only_A(:,:,1:K); % using MUEs old
%% The channel generation is done HERE
% Now is to find the outter precoder such that
%  ctranspose of U * sum of Theta of sue = 0
pathloss_sue = zeros(S,1);
pathloss_sue_dB = zeros(S,1);
% Correlation matrix
for s=1:S
    pathloss_sue_dB(s) = PathLoss(norm(MeNB - HUEnode(:,s)));
    pathloss_sue(s) = 10^(-pathloss_sue_dB(s)/10);    
end
Theta_Sue = zeros(N,N,S);
Theta_Sue_Sum = zeros(N,N); % For Massive MIMO + Wireless
% Theta_Sue_Sum_A = zeros(N,N); % For Massive MIMO + Wired
 for s=1:S
    Theta_Sue(:,:,s) =    US_Var(K+s) * pathloss_sue(s)  * Theta_S(:,:,s) / trace(Theta_S(:,:,s)); 
%     Theta_Sue_A(:,:,s) =    US_Var_A(M+S+s) * pathloss_sue(s) * Theta_S(:,:,s) / trace(Theta_S(:,:,s)); 
    Theta_Sue_Sum = Theta_Sue_Sum + Theta_Sue(:,:,s);
%     Theta_Sue_Sum_A = Theta_Sue_Sum_A + Theta_Sue_A(:,:,s);
 end
  for s=1:S
     HC(:,s) = sqrt(Theta_Sue(:,:,s)) *  HC(:,s);
 end
% Find the outer precoder U
Uprecoder = null(ctranspose(Theta_Sue_Sum));
% Uprecoder_A = null(ctranspose(Theta_Sue_Sum_A));
N_Sc = 0; % If there is a downlink for SUE
for s = K+1:K+S
    if US_Var(s) > 0;
        N_Sc = N_Sc + 1;
    end
end
N_Sc_A = 0;
for s = K+1:K+S
    if US_Var_A(s) > 0;
        N_Sc_A = N_Sc_A + 1;
    end
end
%% Some matrices to calculate the DE of SINR
THETA_tilde_Wireless = zeros(N,N,K);
THETA_tilde_Wired = zeros(N,N,M+S); % # MUEs = # SCs
if N_Sc > 0 % To avoid there is no SC to serve
    for k = 1:K
        THETA_tilde_Wireless(:,:,k) = Uprecoder * ctranspose(Uprecoder) * THETA_tilde_Only(:,:,k)* Uprecoder * ctranspose(Uprecoder);%/trace(Uprecoder * ctranspose(Uprecoder) * THETA(:,:,k)* Uprecoder * ctranspose(Uprecoder)),rounding);
    end
else
        THETA_tilde_Wireless = THETA_tilde_Only;
end
THETA_tilde_Wired = THETA_MUE_Ad;
% if N_Sc_A > 0
%     for m = 1:M+S
%         THETA_tilde_Wired(:,:,m) = Uprecoder_A * ctranspose(Uprecoder_A) * THETA_MUE_Ad(:,:,m)* Uprecoder_A * ctranspose(Uprecoder_A);
%     end
% else
%      THETA_tilde_Wired = THETA_MUE_Ad;
% end
%% 20150506: We here finish the codes to generate the channel propagation, to compute the channel correlation matrix, the outer precoder

%% From here we are going to calculate the cross-tier interference from the SC to the MUE, SC: Have to ignore the thermal noise
IF_S2m = zeros(S,M);
IF_S2s = zeros(S,S);
pathloss_S2m = zeros(S,M);
pathloss_S2m_dB = zeros(S,M);

%% For the MUE
for m = 1:M
    for s=1:S
        ChanG  =   sqrt(1/2)* ( randn + 1i*randn); %Theta_S2u(s,s)  * sqrt(1/2)* ( randn + 1i*randn)
        pathloss_S2m_dB(s,m) = PathLoss(norm(MUEnode(:,m) - HeNBnode(:,s)));
        pathloss_S2m(s,m) = ChanG * 10^(- pathloss_S2m_dB(s,m)/10);
    end
end


%%  For the SC 
pathloss_dB_S2s = zeros(S,S);
pathloss_dB_S2s_dB = zeros(S,S);
for u = 1:S
    for s=1:S
        if s == u
        else
             ChanG  =   sqrt(1/2)* ( randn + 1i*randn); %Theta_S2u(s,s)  * sqrt(1/2)* ( randn + 1i*randn) (N0 -30) +
             pathloss_dB_S2s_dB(u,s) = PathLoss(norm( HeNBnode(:,u) - HeNBnode(:,s))) ;
             pathloss_dB_S2s(u,s) = 10^(-pathloss_dB_S2s_dB(u,s)/10) * ChanG;
        end
    end
end

HH = ctranspose([HM,HS]);
HHH = [HM,HS,HC];
scalefactor = 1; % fucking step to avoid the zero matrix
temp = 0;
for k=1:K
    if US_Var(k) ~=0
%         temp = min(min(min(abs(HH(k,:)))));
        temp = 10^(5)* max(max(max(abs(HH(k,:)))));
    end
    scalefactor= min(scalefactor,temp);
end

% SF1 = (max(max(abs(pathloss_dB_S2s))));
% SF2 = (max(max(abs(pathloss_S2m))));
% NN0 = 10^(-134.5/10);
  
%% SF1=max(SF1,SF2);
%scalefactor = max(scalefactor,SF1);

for m=1:M
    for s=1:S
         IF_S2m(s,m)  =  P_s * norm(pathloss_S2m(s,m)/scalefactor);% - 10*log10(abs(Theta_S2m(s,m))) + 20 * log10(abs(ChanG)); % + 20 * log10(abs(ChanG))
    end
end
for u=1:S
    for s=1:S
         IF_S2s(u,s)  =  P_s * norm(pathloss_dB_S2s(u,s)/scalefactor);
    end
end
% This is interference function from SC to MUEs and other SCs. The format
% is based on the US policy
% scalefactor_Sc = abs(min(min(pathloss_dB_S2m)));
IF_Fun = [ zeros(M,K)
            IF_S2m IF_S2s];
        
I = eye(size(Uprecoder,2),size(Uprecoder,2));
I1 = eye(size(HH,2),size(HH,2));
% HH = HH/scalefactor;
%V = Uprecoder * Tprecoder ;
if N_Sc > 0
    V = Uprecoder * (ctranspose(Uprecoder) * ctranspose(HH) * HH  * Uprecoder / scalefactor^2 +  N * 1 /(100) * I  )^(-1)*ctranspose(Uprecoder) * ctranspose(HH)/ scalefactor^2 ;
else
    V =  (ctranspose(HH) * HH  / scalefactor^2 +  N * 1 /(100) * I1  )^(-1)* ctranspose(HH)/ scalefactor^2 ;
end

%% From here we are going to calculate the SINR from the SC to the correponding SUE
% Theta_S2u = zeros(S,S); % Correlation matrix
SINR_S2u = zeros(1,S);
pathloss_S2u = zeros(S,S);
pathloss_S2u_dB = zeros(S,S);
ChanG_u = zeros(S,S);
for s=1:S
     for u=1:S
        pathloss_S2u_dB(u,s) =  PathLoss(norm(HUEnode(:,u) - HeNBnode(:,s)));
        pathloss_S2u(u,s)= 10^(-pathloss_S2u_dB(u,s)/10) * abs(sqrt(N_s) * ( 0*randn + 1i*1/N_s*randn)).^2;
        ChanG_u(u,s)  =   pathloss_S2u(u,s); 
     end
end
SF3 = (max(max(abs(pathloss_S2u))));


% For open access, we need to calculate the data rate from SC to MUE
SINR_S2M = zeros(M,S);
ChanG_S2M = zeros(M,S);
pathloss_S2M = zeros(M,S);
pathloss_S2M_dB = zeros(M,S);

for m=1:M
    for s=1:S
        pathloss_S2M_dB(m,s) =  PathLoss(norm(MUEnode(:,m) - HeNBnode(:,s)));
        pathloss_S2M(m,s) = 10^(-pathloss_S2M_dB(m,s)/10) ...
            * abs(sqrt(N_s) * (sqrt(1-epsilon^2) * (0/N_s*randn + 1i*1/N_s*randn) + epsilon*(0/N_s*randn + 1i*1/N_s*randn))).^2;
        ChanG_S2M(m,s)  =  pathloss_S2M(m,s);
            
    end
end

SFS2M = max((max(max(abs(pathloss_S2M)))));

for s=1:S
        SINR_S2u(s) = 10^(Tx_Gain/10) * P_s * N_s / N_au * norm( ChanG_u(s,s) )/SF3/(1 +  Epsilon); % + 20 * log10(abs(ChanG)) + 20 * log10(abs(ChanG(s))) In dB, then calculate the rate easily
end

for m =1:M
    for s=1:S
        SINR_S2M(m,s) = 10^(Tx_Gain/10) * P_s * N_s / N_au * norm( ChanG_S2M(m,s) )/SFS2M/((1 + Epsilon)); % + 20 * log10(abs(ChanG)) + 20 * log10(abs(ChanG(s))) In dB, then calculate the rate easily
    end
end

end % end ChannelGene Function

