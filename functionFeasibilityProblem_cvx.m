function [feasible,p_ue,Wsolution] = functionFeasibilityProblem_cvx(H,D,q,gammavar)
%Solves the feasibility problem with quality-of-service (QoS) constraints for uplink.
%this function using the cvx to solve the virtual downlink problem.
%The power minimization under QoS requirements is
%
%minimize   betavar
%subject to SINR_k >= gammavar(k) for all users k,
%           Power constraints scaled by betavar.
%
%If this optimization problem is feasible and betavar<=1, then the
%feasibility problem with QoS constraints is also feasible.
%
%This optimization problem is convex. The computational complexity is
%therefore polynomial in the number of users, antennas, and power
%constraints. The implementation can, at least, handle 30 users, 50
%antennas, and 50 power constraints.
%
%This is version 1.1.
%The implementation utilizes and requires CVX: http://cvxr.com/
%
%INPUT:
%H          = Kt x Kt*Nt matrix with row index for receiver and column
%             index transmit antennas
%D          = Kt*Nt x Kt*Nt x Kt diagonal matrix. Element (j,j,k) is one if
%             j:th transmit antenna can transmit to user k and zero otherwise
%q          = Limits of the L power constraints
%gammavar   = Kt x 1 vector with SINR constraints for all users.
%
%OUTPUT:
%feasible  = This variable is feasible=true if the feasibility problem is
%            feasible. Otherwise we have feasible=false.
%Wsolution = Kt*Nt x Kt matrix with beamforming achieved by the power
%            minimization problem. This matrix is empty if this problem is
%            infeasible.


K = size(H,1); %Number of users
N = size(H,2); %Number of transmit antennas (in total)
L = K; %Number of links

%Channels of the signal intended for user i when it reaches user k
hkD = zeros(K,N);
for i = 1:K
    hkD(i,:) = H(i,:)*D(:,:,i);
end

%Solve the power minimization under QoS requirements problem using CVX
cvx_begin
cvx_quiet(true); % This suppresses screen output from the solver

variable Wv(N,K) complex
variable t(K) nonnegative;
dual variables y{K}

minimize sum(t) %Minimize the power indirectly by scaling power constraints

subject to


%SINR constraints (K constraints)
for k = 1:K  
    
    norm([Wv(:,k)]) <= sqrt(t(k));
    
    imag(hkD(k,:)*Wv(:,k)) == 0; %Useful link is assumed to be real-valued
    
    %SOCP formulation for the SINR constraint of user k
    y{k}: real(hkD(k,:)*Wv(:,k)) >= sqrt(gammavar(k))*norm([1 hkD(k,:)*Wv(:,[1:k-1 k+1:K])]);  
end

% for downlink
% %Power constraints (L constraints) scaled by the variable betavar
% for l = 1:L
%     norm( Wv(:,l)*Qsqrt,'fro') <= betavar*sqrt(q);
% end
% 
% betavar >= 0; %Power constraints must be positive

%for v-downlink
cvx_end
for l=1:L
    P(l,1) = norm( Wv(:,l),'fro')^2;
end
betavar=max(  P./q );

%Analyze result and prepare the output variables.
if isempty(strfind(cvx_status,'Solved')) %Both power minimization problem and feasibility problem are infeasible.
    feasible = false;
    p_ue = q;
    Wsolution = [];
elseif betavar>1 %Only power minimization problem is feasible.
    feasible = false;
    p_ue = q;
    Wsolution = Wv;
else %Both power minimization problem and feasibility problem are feasible.
    feasible = true;
    p_ue = P ;
    Wsolution = Wv;
end
