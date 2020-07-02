%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% WSINDy: function for recovering weights using SINDy
%%%%%%%%%%%% 
%%%%%%%%%%%% ode_num selects an ODE system from the list ode_names
%%%%%%%%%%%% tol_ode sets the tolerance (abs and rel) of ode45 
%%%%%%%%%%%% noise_ratio sets the signal-to-noise ratio (L2 sense)
%%%%%%%%%%%% 
%%%%%%%%%%%% Copyright 2020, All Rights Reserved
%%%%%%%%%%%% Code by Daniel A. Messenger
%%%%%%%%%%%% For Paper, "Weak SINDy: Galerkin-based Data-Driven Model
%%%%%%%%%%%% Selection"
%%%%%%%%%%%% by D. A. Messenger and D. M. Bortz

function w_sparse_sindy = standard_sindy(t,xobs,Theta_0,M_diag, useFD,n,lambda,gamma)

if useFD==0
    useFD=1;
    L = length(t)-useFD*2;
    C = fdcoeffF(1, t(useFD+1), t(1:useFD*2+1));
    D = spdiags(repmat(C(:,end)',L,1), 0:useFD*2 ,L,length(t));
    naive_diff = D*xobs;
    dt = mean(diff(t));
    dxobs_0 = 0*naive_diff;
    for j=1:n
        dxobs_0(:,j) = TVRegDiff( xobs(2:end-2,j), 10, dt, naive_diff(:,j), [], [], dt, 0, 0 );
    end
    useFD=1;
else
    L = length(t)-useFD*2;
    C = fdcoeffF(1, t(useFD+1), t(1:useFD*2+1));
    D = spdiags(repmat(C(:,end)',L,1), 0:useFD*2 ,L,length(t));
    dxobs_0 = D*xobs;
end

if ~isempty(M_diag)
    w_sparse_sindy = sparsifyDynamics(Theta_0(useFD+1:end-useFD,:).*(1./M_diag'),dxobs_0,lambda,n,gamma);
    w_sparse_sindy = (1./M_diag).*w_sparse_sindy;
else
    w_sparse_sindy = sparsifyDynamics(Theta_0(useFD+1:end-useFD,:),dxobs_0,lambda,n,gamma);
end    
end
