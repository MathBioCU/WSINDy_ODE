    
function [Xi,its] = sparsifyDynamics(Theta,dXdt,lambda,gamma)
    n = size(dXdt,2);
    nn =size(Theta,2);
    if  gamma ~= 0
        Theta = [Theta;gamma*eye(nn)];
        dXdt = [dXdt;zeros(nn,n)];
    end
    Xi = Theta \ dXdt;  % initial guess: Least-squares
    smallinds = 0*Xi;
    for j=1:nn
        smallinds_new = (abs(Xi)<lambda);
        if all(smallinds_new(:)==smallinds(:))
            its = j;
            return
        else
            smallinds = smallinds_new;
            Xi(smallinds)=0;
            for ind = 1:n        
                biginds = ~smallinds(:,ind);
                Xi(biginds,ind) = Theta(:,biginds)\dXdt(:,ind);
            end
        end
    end
    its = nn;
end
