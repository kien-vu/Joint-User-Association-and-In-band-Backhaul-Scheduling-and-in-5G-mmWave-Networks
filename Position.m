function [ MUEnode,  MUEnode_Ad, HeNBnode, HUEnode] = Position( M, S)
%  This function is used to generate the eNodeB and users position.
%  The Covariance matrix also done here
%  Node Position Generation Square Area
global Axis; % in order change the network size 
% The maximum transmit power should change accordingly to the area. Based
% on the SNR of cell-edge user
global P_b0;
global Tx;  % Total transit power at the SC
%%
    % 20 MHz -101 dBm
    % 100 MHz -94.1 dBm
    % 1000 MHz -84 dBm
%%
X_Axis = Axis; % the size of network area
Y_Axis = Axis;
global SCperRow % Number of SCs per each row or column
HeNBnode = rand(2,S);
MUEnode  = rand(2,M);
MUEnode_Ad = rand(2,S);
HUEnode = rand(2,S);
Original = [ X_Axis/2; X_Axis/2];
global Bound1 % To ensure the MUEs are located faw way MBS with distance Bound1, and HUEs are located inside the coverage of SCs
global Bound2 % To ensure the SCs are located faw way MBS with distance Bound2
Bound1 = 35; 
Bound2 = 7.5; % To ensure the MUEs are located faw way SCs with distance Bound2

% Here we locate the SCs into a grid
  for x = 1:SCperRow
      for y = 1:SCperRow
          HeNBnode(1,(x-1)*SCperRow + y) = X_Axis/SCperRow *1/2 +(x-1)*X_Axis/SCperRow;
          HeNBnode(2,(x-1)*SCperRow + y) = X_Axis/SCperRow *1/2 +(y-1)*X_Axis/SCperRow;
            while(1)
                if norm(HeNBnode(:,(x-1)*SCperRow + y) - Original) <= Bound1
                            HeNBnode(1,(x-1)*SCperRow + y) =   HeNBnode(1,(x-1)*SCperRow + y) + 1 * (-5+10*randn);
                            HeNBnode(2,(x-1)*SCperRow + y) =   HeNBnode(2,(x-1)*SCperRow + y) + 1 * (-5+10*randn);
                else 
                    break;
                end
            end
      end
  end
% Here small cell users are located surround the SCs
  for s = 1:S
        HUEnode(1,s) =  HeNBnode(1,s) + 1 * (-5+10*rand);
        HUEnode(2,s) =  HeNBnode(2,s) + 1 * (-5+10*rand); 
             while(1)
                if norm(HUEnode(:,s) - HeNBnode(:,s)) > Bound1
                    HUEnode(1,s) =  HeNBnode(1,s) + 1 * (-5+10*rand);
                    HUEnode(2,s) =  HeNBnode(2,s) + 1 * (-5+10*rand); 
                else 
                    break;
                end
             end
            while(1)
                if norm(HUEnode(:,s) - HeNBnode(:,s)) < 2
                    HUEnode(1,s) =  HeNBnode(1,s) + 1 * (-5+10*rand);
                    HUEnode(2,s) =  HeNBnode(2,s) + 1 * (-5+10*rand); 
                else 
                    break;
                end
            end
  end
% Here for MUEs
   for m = 1:M
            MUEnode(1,m) =   X_Axis * rand;
            MUEnode(2,m) =   Y_Axis * rand;
            while(1)
                if norm(MUEnode(:,m) - Original) <= Bound1
                            MUEnode(1,m) =   X_Axis * rand;
                            MUEnode(2,m) =   Y_Axis * rand;
                else 
                    break;
                end
            end
            for s=1:S
                while(1)
                    if norm(MUEnode(:,m) - HeNBnode(:,s)) <= Bound2
                                MUEnode(1,m) =   X_Axis * rand;
                                MUEnode(2,m) =   Y_Axis * rand;
                    else 
                        break;
                    end
                end
            end
   end  
     for s=1:S
        MUEnode_Ad(1,s) = X_Axis * rand;
        MUEnode_Ad(2,s) = X_Axis * rand;
            while(1)
                if norm(MUEnode_Ad(:,s) - Original) <= Bound1
                            MUEnode_Ad(1,s) =   X_Axis * rand;
                            MUEnode_Ad(2,s) =   Y_Axis * rand;
                else 
                    break;
                end
            end
     end   
%  figure;
%  plot(Original(1,:),Original(2,:),'r-o'); hold on
%  plot(MUEnode(1,:),MUEnode(2,:),'ks'); hold on 
%  plot(HeNBnode(1,:),HeNBnode(2,:),'b*'); hold on 
%  plot(HUEnode(1,:),HUEnode(2,:),'g+'); hold on 



end