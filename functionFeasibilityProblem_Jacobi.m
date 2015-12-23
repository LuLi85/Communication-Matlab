function [feasible,p_ue,Wsolution] = functionFeasibilityProblem_Jacobi(H,Gamm,pmax,NC)
%this function solves the new feasibility problem including BSs assignment
%with quality-of-service (QoS) constraints
%INPUT:
%H          = K x L*Nt matrix with row index for receiver and column
%             index transmit antennas, L is the total BSs
%A          = K x L . Element (k,l) is one if
%             the l-th transmitter are assigned to user k and zero otherwise
%Gamm       = K x 1  vector with Gamm(k) the SINR threshold
%pmax       = Limits of the Kr power constraints P_max
%
%NC         = scalar the constant number of BSs

%OUTPUT:
%feasible  = This variable is feasible=true if the feasibility problem is
%            feasible. Otherwise we have feasible=false.
%Wsolution = L*Nt x Kr matrix with beamforming achieved by the power
%            minimization problem. This matrix is empty if this problem is
%            infeasible.
global LBS;
K = size(H,1); %Number of users
N = size(H,2); %Number of transmit antennas (in total)
Nt = N/LBS; %Number of transmit antennas of each transmitter
global mode
if (mode==0)
    L=NC; %the number of links when NC active BSs
elseif (mode==1) %1:the user-centric BSs assignment
    L=K; %the number of equaivalent links
end

%init
w=ones(LBS*Nt,K)/sqrt(LBS*Nt); % init normailized w randomly
P=pmax.*ones(K,1); % init transmit Power
tol = 0.0001; %error tolorence
no = ones(K,1);%noise variance

%generate the search space of BSs assignment for NC=2
BS_selection;

%start traversing
for s = 1: size(A_F,3)
    A = A_F(:,:,s);
    
    A=kron(A,ones(1,Nt));
    %obtain the effective channel a
    G = H';
    %the interference from user j to i
    for i=1:1:K
        for j=1:1:K
            a{i,j} = G(:,j).*A(i,:)';
        end
    end
    
    %algorithm starts here
    iterations=1;
    Err=1; %some initial error value
    
    while max(Err)>tol & iterations<30   % I choose maximum erro to be a divergence criteria
        %if number of transmit antennas is greater than number of user it must be feasible
        %update w(n+1) for each user
        for k=1:K
            temp=zeros(LBS*Nt,LBS*Nt);
            for j=1:K
                if (j~=k)
                    temp = P(j)*a{k,j}*a{k,j}'+ temp;
                end
            end
            temp=temp+eye(LBS*Nt,LBS*Nt)*no(i); %noise is one
            % temp=pinv(temp);
            w(:,k)= pinv(temp)*a{k,k};
        end
        
        %update G(w(n+1))
        for k=1:K
            for j=1:K
                if j~=k
                    G_w(k,j)=abs(w(:,k)'*a{k,j})^2/abs(w(:,k)'*a{k,k})^2;
                end
            end
        end
        
        %update P(n+1) for each user
        for k=1:L
            temp=0;
            for j=1:L
                if(j~=k)
                    temp = temp+ G_w(k,j)*P(j);
                end
            end
            P(k) = Gamm(k)*temp + Gamm(k)*  no(k)* abs(w(:,k)'*w(:,k))/ (abs(w(:,k)'*a{k,k})^2);
            % P(k) = min([P(k), pmax(k)]);
            % P(k) = max([0, P(k)]);
        end
        
        iterations=iterations+1;
        
        %     for i=1:K
        %         SINR(i,iterations)=P(i)/(G_w(i,:)*P+no(i)* (w(:,i)'*w(:,i))/(abs(w(:,i).'*a{i})^2) );
        %     end
        %     Err=abs(Gamm- SINR(:,iterations)) %error
        
    end
    
    %check if P satisfies 0< P < Pmax
    if  ( sum( P > pmax)>=1 | sum( P < 0)>=1 )
        feasible = false;
        Wsolution = [];
        p_ue = pmax;
    else
        %change the feasibility state and break off
        feasible = true;
        Wsolution = w;
        p_ue = P;
        break;
    end
    
end
