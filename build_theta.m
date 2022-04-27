%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% WSINDy: function for building data matrix Theta_0
%%%%%%%%%%%% 
%%%%%%%%%%%% Copyright 2020, All Rights Reserved
%%%%%%%%%%%% Code by Daniel A. Messenger
%%%%%%%%%%%% For Paper, "Weak SINDy: Galerkin-based Data-Driven Model
%%%%%%%%%%%% Selection"
%%%%%%%%%%%% by D. A. Messenger and D. M. Bortz

function [Theta_0, tags, true_nz_weights, M_diag, lambda] = build_theta(xobs,weights,common_params,n)

%%%%%%%%%%%%%%%%% get params

polys = common_params{1}; 
trigs = common_params{2}; 
lambda_mult = common_params{3};
scale_Theta = common_params{4};

%%%%%%%%%%%%%%%%% build Theta_0 & true_nz_weights

[Theta_0,tags] = poolDataGen(xobs,polys,trigs);
true_nz_weights = get_true_weights(weights,tags,n);

if scale_Theta > 0
    M_diag = vecnorm(Theta_0,scale_Theta,1)'; 
    true_weights_scale = M_diag.*true_nz_weights; 
    lambda = min(abs(true_weights_scale(true_nz_weights~=0)))/lambda_mult;
else
    M_diag = [];
    lambda = min(abs(true_nz_weights(true_nz_weights~=0)))/lambda_mult; 
end

end

function [Theta_0,tags] = poolDataGen(xobs,polys,trigs)

n = size(xobs,2);
ind = 0;

tags = [];
for p = 1:length(polys)
    monom_powers = partitionNk(polys(p),n);
    num_monoms = size(monom_powers,1);
    for j = 1:num_monoms
        Theta_0(:,ind+1) = prod(xobs.^(monom_powers(j,:)),2);
        ind = ind+1;
    end
    tags = [tags;monom_powers];
end
for k=1:length(trigs)
    trig_inds = [-trigs(k)*1i*eye(n);trigs(k)*1i*eye(n)];
    Theta_0 = [Theta_0 sin(trigs(k)*xobs) cos(trigs(k)*xobs)];
    tags = [tags; trig_inds];
end
end

function total_parts = partitionNk(N,k)
    if N>1
%         A = partitions(N,0:N,N,k);
%         num_unique_parts = size(A,1);
%         unique_parts = zeros(num_unique_parts,k);
%         for j=1:num_unique_parts
%             temp_part = [];
%             for m= 1:N+1
%                 if A(j,m)>0
%                     for l = 1:A(j,m)
%                         temp_part = [temp_part m-1];
%                     end
%                 end
%             end
%             unique_parts(j,:) = temp_part;
%         end
%         total_parts = [];
%         for j=1:num_unique_parts
%             total_parts = [total_parts;unique(perms(unique_parts(j,:)),'rows')];
%         end
        total_parts = partitions(N,ones(1,k));
    elseif N==1
        total_parts = eye(k);
    elseif N==0
        total_parts = zeros(1,k);
    end
end

function true_nz_weights = get_true_weights(weights,tags,n)
    true_nz_weights = zeros(size(tags,1),n);
    for i = 1:length(weights)
        weights_i = weights{i};
        [l1,l2] = size(weights_i);
        for j = 1:l1
            true_nz_weights(all(weights_i(j,1:l2-1) == tags,2),i) = weights_i(j,l2);
        end
    end
end