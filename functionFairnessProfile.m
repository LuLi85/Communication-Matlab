function [finalInterval,WBestBeamforming,nbrOfEvaluations] = functionFairnessProfile(H,D,q,delta,lowerPoint,upperPoint,feasibilityMode)
%This is a generalization of max-min fairness optimization. The FPO
%problem searches on a line between LOWERPOINT and UPPERPOINT and finds
%the intersection between the line and the Pareto boundary of the rate
%region. The search is based on bisection and solving a power minimization

%under QoS requirements;
%
%The FPO problem is equivalent to finding (g_1,...,) that solves
%
%maximize min_k ( g_k - lowerPoint(k) ) / ( upperPoint(k) - lowerPoint(k) )
%
%subject to     g_k >= lowerPoint(k) for all users k,
%               Power constraints.
%
%This problem is quasi-convex, meaning that it can be solved as a short
%sequence of convex optimization problem. The computational complexity is
%therefore polynomial in the number of users, antennas, and power
%constraints.
%This is version uplink.
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%original article listed above.
%
%The implementation utilizes and requires CVX: http://cvxr.com/
%
%INPUT:
%H           = NC x Kt*Nt matrix with row index for receiver and column
%              index transmit antennas
%D           = Kt*Nt x Kt*Nt x Kr diagonal matrix. Element (j,j,k) is one if
%              j:th transmit antenna can transmit to user k and zero otherwise
%q           = Limits of the L power constraints
%delta       = Accuracy of the final solution. The algorithm terminates when
%                 norm(upperPoint - lowerPoint) <= delta
%lowerPoint  = Start point of the line (must be inside the rate region)
%upperPoint  = End point of the line (must be outside of the rate region)
%feasiblityMode = (Optional) Consider different feasibility checking
%algorithm

%OUTPUT:
%finalInterval    = Kr x 2 matrix with lowerPoint and upperPoint at
%                   termination
%WBestBeamforming = Kt*Nt x Kr matrix with beamforming that achieves the
%                   lower point in the final interval
%nbrOfEvaluations = Number of times that the convex subproblem (power
%                   minimization under QoS requirements) is solved

if nargin < 8
    specialMode = 0;
end

global LBS;
K = size(H,1); %Number of users
NC = 2; %the constant number of BSs
Nt = size(H,2)/LBS;

%Pre-allocation of matrix for storing optimal beamforming
WBestBeamforming = [];

%Count the number of feasibility problem solved
nbrOfEvaluations = 0;

%%Part 1: Solve the problem by bisection.

%Solve the problem by bisection - iterate until different between
%current lower and upper point
while norm(upperPoint - lowerPoint) > delta
    
    candidatePoint = (lowerPoint+upperPoint)/2; %Compute the midpoint at the line
    
    gammavar = 2.^(candidatePoint)-1; %Transform midpoint into SINR requirements
    
    %Check the feasibility at the midpoint by solving a feasibility problem
    if (feasibilityMode == 1)
        [feasible,p_ue,Wcandidate] = functionFeasibilityProblem_cvx(H,D,q,gammavar);
        % under dynamic BSs assignment
    elseif (feasibilityMode == 2)
        [feasible,p_ue,Wcandidate] = functionFeasibilityProblem_Jacobi(H,gammavar,q,NC);
    end
    
    %If the problem was feasible, then replace lowerPoint with
    %candidatePoint and store W as current best solution.
    if feasible
        lowerPoint = candidatePoint;
        WBestBeamforming = Wcandidate;
        PU = p_ue;
    else
        %If the problem was not feasible,then replace upperPoint with candidatePoint
        upperPoint = candidatePoint;
    end
    
    %Increase the number of function evaluations
    nbrOfEvaluations = nbrOfEvaluations+1;
end




%%Part 2: Prepare the achieved solution for output
%If the midpoints analyzed by the algorithm have never been feasible,
%then obtain a feasible beamforming solution using the lowerPoint. This
%happens when delta is too large or when the optimal point is very
%close to lowerPoint.
if isempty(WBestBeamforming)
    gammavar = 2.^(lowerPoint)-1;
    
    %Check the feasibility at the midpoint by solving a feasibility problem
    if (feasibilityMode == 1)
        [feasible,p_ue,Wcandidate] = functionFeasibilityProblem_cvx(H,D,q,gammavar);
        % under dynamic BSs assignment
    elseif (feasibilityMode == 2)
        [feasible,p_ue,Wcandidate] = functionFeasibilityProblem_Jacobi(H,gammavar,q,NC);
    end
    
    if feasible
        WBestBeamforming = Wcandidate;
    else
        %The algorithm requires that the start point is inside of the
        %rate region, which is not the case if we end up here.
        error('Fairness-profile optimization problem is infeasible');
    end
end

%Prepare for output
%the uplink channel
G=H';
save WBestBeamforming;
%calculate the rate
for k=1:K
    %the effective signal power
    cg(k) =  abs(WBestBeamforming(:,k)'*G(:,k))^2*PU(k);
    %the receive noise power
    ig(k) = WBestBeamforming(:,k)'*eye(LBS*Nt)*WBestBeamforming(:,k);
    %the interference power
    for j=1:K
        if (j~=k)
            ig(k) = ig(k)+abs(WBestBeamforming(:,k)'*G(:,j))^2*PU(j);
        end
    end
    %the output rates
    rates(k,1) = log2(1+cg(k)/ig(k));
end
%   Compute the rates that are actually achieved
%     channelGains = abs(WBestBeamforming*H(k,:)).^2;
%     signalGains = diag(channelGains);
%     interferenceGains = sum(channelGains,2)-signalGains;
%     rates = log2(1+signalGains./(1+interferenceGains));

%Store the final interval between lower and upper point
if sum(rates>lowerPoint) == K
    finalInterval = [rates upperPoint];
else
    finalInterval = [lowerPoint upperPoint];
end

