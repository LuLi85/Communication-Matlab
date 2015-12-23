%new parameters:
%NC: the selected BSs for serving
%mode: network-centric / user-centric selection

%the global BSs assignment matrix
A = ones(K,LBS);
%determine A
if(NC==2)
    index=[];tmp=zeros(K,1);
    % select NC BSs to be associated
    for s1 =1: LBS
        for s2 =1:LBS
            tmp(1) = s1;% the first BS
            tmp(2) = s2;% the second BS
            % if they are not the same one
            if (tmp(1) ~= tmp(2) & tmp(1)<tmp(2) )
                %save into the set of selected BSs
                index = [index, tmp];
            end
        end
    end
    
    %network-centric BSs assignment
    if (mode==0)
        for i=1:size(index,2)
            A=zeros(K,LBS);%reset A
            for s=1:NC
                %every selected BS can serve two users
                A(:,index(s,i)) =ones(K,1);
            end
            %A_F is A-matrices 
            A_F(:,:,i) = A;
        end
    end
    %user-centric BSs assignment
    if(mode==1)
        %for user 1
        for i1=1:size(index,2)
            %for user 2
            for i2=1:size(index,2)
                A=zeros(K,LBS);%reset A
                for s=1:NC
                    A(1,index(s,i1)) =1;
                    A(2,index(s,i2)) =1;
                end
                %A_F is A-matrices
                A_F(:,:,(i1-1)*size(index,2)+i2) = A;
            end
        end
    end
    
end
