%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% WSINDy: script for recoverying various ODE systems
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


%%% generate clean data 

ode_num = 6;                        % select ODE system from the list ode_names
tol_ode = 1e-10;                    % set tolerance (abs and rel) of ode45
noise_ratio = 0.1;                  % set signal-to-noise ratio (L2 sense)
rng('shuffle');rng_seed =rng; rng_seed =rng_seed.Seed;

ode_names = {'Linear','Logistic_Growth','Van_der_Pol','Duffing','Lotka_Volterra','Lorenz'};
ode_name = ode_names{ode_num};

if strcmp(ode_name, 'Linear')
    ode_params = {[[-0.1 2];[-2 -0.1]]};x0 = [3;0]; tspan = 0:0.01:15;
    
    % set params common to WSINDy and SINDy 
    polys = 0:5;                    % monomial powers to include in library
    trigs = [];                     % sine / cosine frequencies to include in library
    gamma = 10^-Inf;                % Tikhonoff regularization parameter
    lambda_mult = 4;                % sets sparsity knob lambda = min(true_weights)/lambda_mult; lambda_mult = Inf => lambda = 0
    scale_Theta = 0;                % toggle normalize columns of Theta_0 by norm = scale_Theta_0

    % set WSINDy params            
    tau = 1;                        % toggle adaptive grid - convex comb parameter between uniform (=0) and adaptive (=1)
    K = 126;                        % num basis fcns
    p = 2; s = 16;                  % AG weak-deriv params -- poly degree, support size
    r_whm = 30;                     % width at half-max in num of timepoints (s.t. phi(r_whm*dt) = 1/2) 
    useGLS = 1;                     % toggle generalized least squares

elseif strcmp(ode_name, 'Logistic_Growth')
    ode_params = {2}; x0 = 0.01; tspan = 0:0.005:10;

    % set common params
    polys = 0:5;                    
    trigs = [];                    
    gamma = 10^-2.0;                
    lambda_mult = 4;                
    scale_Theta = 2;               

    % set WSINDy params             
    tau = 1;                        
    K = 100;                        
    p = 2; s = 16;                  
    r_whm = 30;                     
    useGLS = 1;                     
    
elseif strcmp(ode_name, 'Van_der_Pol')
    ode_params = {4}; x0 = [0;1]; tspan = 0:0.01:30;    

    % set common params
    polys = 0:5;                   
    trigs = [];                     
    gamma = 10^-Inf;
    lambda_mult = 4;
    scale_Theta = 0;

    % set WSINDy params
    tau = 1;
    K = 126;
    p = 2; s = 16;
    r_whm = 30; 
    useGLS = 1;

elseif strcmp(ode_name, 'Duffing')
    mu = 0.2;ode_params = {mu, mu^2/4*5,1}; x0 = [0;2]; tspan = 0:0.01:30;       

    % set common params
    polys = 0:5;
    trigs = [];
    gamma = 10^-Inf;
    lambda_mult = 2;
    scale_Theta = 0;
    
    % set WSINDy params
    tau = 1;
    K = 126;
    p = 2; s = 16;
    r_whm = 30;
    useGLS = 1;

elseif strcmp(ode_name, 'Lotka_Volterra')
    alpha=2/3;beta = 4/3; ode_params = {alpha,beta,1,1}; x0 = [10;10]; tspan = 0:0.02:200;

    % set common params
    polys = 0:5;
    trigs = [];
    gamma = 0;
    lambda_mult = 4;
    scale_Theta = 0;

    % set WSINDy params
    tau = 1;
    K = 126;
    p = 2; s = 16;
    r_whm = 30;
    useGLS = 1;

elseif strcmp(ode_name, 'Lorenz')
    ode_params = {10, 8/3, 28}; tspan = .001:.001:10;
    x0 = [-8 7 27]';
    %x0 = [rand(2,1)*30-15;rand*30+10]

    % set common params
    polys = 1:5;
    trigs = [];
    gamma = 10^-Inf;
    lambda_mult = 4;
    scale_Theta = 0;
    
    % set WSINDy params
    tau = 1;
    K = 224;
    p = 2; s = 16;
    r_whm = 30;
    useGLS = 1;
end

n = length(x0); 
[weights,t,x,rhs] = sim_ode(x0,tspan,tol_ode,ode_name,ode_params);

%%% add noise

signal_power = rms(x(:));
sigma = noise_ratio*signal_power;
noise = normrnd(0,sigma,size(x));
xobs = x + noise;
tobs = t;

%%% recover dynamics

tic;
common_params = {polys,trigs,lambda_mult,scale_Theta,gamma};
wsindy_params = {s, K, p, tau};

[Theta_0, tags, true_nz_weights, M_diag, lambda] = build_theta(xobs,weights,common_params,n);

w_sparse = zeros(size(Theta_0,2),n);
mats = cell(n,1);
ps_all = [];
ts_grids = cell(n,1);
RTs = cell(n,1);
Ys = cell(n,1);
Gs = cell(n,1);
bs = cell(n,1);

for i=1:n
    [Y,grid_i] = adaptive_grid(tobs,xobs(:,i),wsindy_params);
    [V,Vp,ab_grid,ps] = VVp_build_adaptive_whm(tobs,grid_i, r_whm, {0,inf,0});  %{pow,nrm,ord}. ord=0, ||phi||, ord=1, ||phi'|| 
    ps_all = [ps_all;ps];
    mats{i} = {V,Vp};
    ts_grids{i} = ab_grid;
    Ys{i} = Y;
    
    if useGLS == 1
        Cov = Vp*Vp'+10^-12*eye(size(V,1));
        [RT,flag] = chol(Cov);        
        RT = RT';
        G = RT \ (V*Theta_0);
        b = RT \ (Vp*xobs(:,i));
    else
        RT = (1./vecnorm(Vp,2,2));
        G = V*Theta_0.*RT;
        b = Vp*xobs(:,i).*RT;
    end        
   RTs{i} = RT;
   Gs{i} = G;
   bs{i} = b;
 
    if scale_Theta > 0
        w_sparse_temp = sparsifyDynamics(G.*(1./M_diag'),b,lambda,1,gamma);
        w_sparse(:,i) = (1./M_diag).*w_sparse_temp;
    else
        w_sparse(:,i) = sparsifyDynamics(G,b,lambda,1,gamma);
    end
end
ET = toc;

%%% get SINDy model

useFD = 0;                      % finite difference differentiation order, if =0, then uses TVdiff
w_sparse_sindy = standard_sindy(t,xobs,Theta_0,M_diag, useFD,n,lambda,gamma);
err_sindy = [norm(w_sparse_sindy-true_nz_weights);norm(w_sparse_sindy(true_nz_weights~=0)-true_nz_weights(true_nz_weights~=0))]/norm(true_nz_weights);        

%%% visualize basis and covariance, display error analysis in command window

figure(3); clf
set(gcf, 'position', [1250 10 700 450])
for d=1:n
subplot(3,n,d) 
plot(tobs,xobs(:,d),'r-',tobs(floor(mean(ts_grids{d},2))),mean(xobs(:,d))*ones(size(ts_grids{d})),'.k')
subplot(3,n,n+d)
plot(tobs,mats{d}{1}')
subplot(3,n,2*n+d)
spy(RTs{d})
end
 
supp_range = [];
for i=1:n    
    supp_range = [supp_range; ts_grids{i}*[-1;1]];
end

clc;
disp(['log10 2norm err (all weights) (WSINDy)=',num2str(log10(norm(true_nz_weights(:)-w_sparse(:))/norm(true_nz_weights(:))))])
disp(['log10 2norm err (all weights) (SINDy)=',num2str(log10(err_sindy(1)))])
disp(['log10 2norm err (true nz weights) (WSINDy)=',num2str(log10(norm(true_nz_weights(true_nz_weights~=0)-w_sparse(true_nz_weights~=0))/norm(true_nz_weights)))])
disp(['log10 2norm err (true nz weights) (SINDy)=',num2str(log10(err_sindy(2)))])
disp(['min/mean/max deg=',num2str([min(ps_all) mean(ps_all) max(ps_all)])])
disp(['min/mean/max supp=',num2str([min(supp_range) mean(supp_range) max(supp_range)])])
disp(['Num Basis elts=',num2str(length(ps_all))])
disp(['Run time=',num2str(ET)])
