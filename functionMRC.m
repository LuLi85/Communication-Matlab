function wMRC = functionMRC(G,D)
%Calculates the maximum ratio combination (MRC) beamforming vectors.
%for the scenario where coordinated BSs are associated to one user.
%
%This is version 1.0.
%
%INPUT:
%G  = LBS*Nt x K matrix with row index for users and column index
%     transmit antennas
%D  = LBS*Nt x LBS*Nt x K diagonal matrix. Element (j,j,k) is one if j:th
%     receive antenna can receive user k and zero otherwise
%
%OUTPUT:
%wMRT = Kt*Nt x Kr matrix with normalized MRC beamforming



%Number of users
K = size(G,2);

%Total number of total antennas
N = size(G,1);

%Pre-allocation of MRC beamforming
wMRC = zeros(size(G'));

%Computation of MRC, based on Definition 3.2
for k = 1:K
    channelvector = (G(:,k))'; %Useful channel
    wMRC(k,:) = channelvector/norm(channelvector); %Normalization of useful channel
end
