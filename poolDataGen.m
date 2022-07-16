
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% WSINDy: function for building data matrix Theta_0
%%%%%%%%%%%%
%%%%%%%%%%%% Copyright 2020, All Rights Reserved
%%%%%%%%%%%% Code by Daniel A. Messenger
%%%%%%%%%%%% For Paper, "Weak SINDy: Galerkin-based Data-Driven Model                                
%%%%%%%%%%%% Selection"
%%%%%%%%%%%% by D. A. Messenger and D. M. Bortz

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
