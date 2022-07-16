%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% WSINDy_PDE: solve regularized LSP with sequential thresh.
%%%%%%%%%%%% 
%%%%%%%%%%%% Copyright 2020, All Rights Reserved
%%%%%%%%%%%% Code by Daniel A. Messenger
%%%%%%%%%%%% For Paper, "Weak SINDy for Partial Differential Equations"
%%%%%%%%%%%% by D. A. Messenger and D. M. Bortz

function [W,lossvals,its] = wsindy_pde_RGLS_seq(lambdas,gamma,Theta_pdx,lhs_ind,M_scale,alpha)

num_eq = length(lhs_ind);
[K,m] = size(Theta_pdx);
    
G = Theta_pdx(:,~ismember(1:m,lhs_ind));
b = zeros(K,num_eq);
for k=1:num_eq
    b(:,k) = Theta_pdx(:,lhs_ind(k));
end
%--------------------------------------------------------------
W_ls = G \ b;
GW_ls = norm(G*W_ls);

proj_cost = [];
overfit_cost = [];
lossvals = [];

W = zeros(m-num_eq,num_eq);
for l=1:length(lambdas)
    lambda = lambdas(l);
    M = [];
    for k=1:num_eq
        if isempty(M_scale)
            [W(:,k),~] = sparsifyDynamics(G, b(:,k), lambda, 1, gamma,[]);
        else
            M = [M M_scale(~ismember(1:m,lhs_ind))/M_scale(lhs_ind(k))];
            [W(:,k),~] = sparsifyDynamics(G, b(:,k), lambda,1, gamma,M(:,end));
            W(:,k)=W(:,k)./M(:,end);
        end
    end
    
    proj_cost = [proj_cost 2*alpha*norm(G*(W-W_ls))/GW_ls];
%    proj_cost = [proj_cost 2*alpha*abs(norm(G*W-b)-norm(G*W_ls-b))/GW_ls];
    overfit_cost = [overfit_cost 2*(1-alpha)*length(find(W~=0))/length(W(:))];
    lossvals = [lossvals proj_cost(end) + overfit_cost(end)];
end

l = find(lossvals == min(lossvals),1);

lambda = lambdas(l);

M = []; its = [];
for k=1:num_eq
    if ~isempty(M_scale)
        M = [M M_scale(~ismember(1:m,lhs_ind))/M_scale(lhs_ind(k))];
        [W(:,k),its(k)] = sparsifyDynamics(G, b(:,k), lambda,1, gamma,M(:,end));
    else
        [W(:,k),its(k)] = sparsifyDynamics(G, b(:,k), lambda, 1, gamma,[]);
    end
end

lossvals = [lossvals;lambdas; [[lossvals(1:l);lambdas(1:l)] zeros(2,length(lambdas)-l)]; proj_cost; overfit_cost];
end
