clc;
% close all; 
% clear all;

%% Generate Data (xobs,tobs,weights)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ode_num = 4;                        % select ODE from list below
tol_ode = 1e-12;                    % set tolerance (abs and rel) of ode45
noise_ratio = 0.1;                    % set ratio of noise to signal (L2 sense), if negative, sets sigma to -noise_ratio
rng('shuffle');                   % comment out to reproduce previous noise 
rng_seed =rng().Seed;
rng(rng_seed);

dt = []; tspan = [0:0.01:30]; ode_params = {0.1, 0.1, 5}; x0 = [0;1]; % change as desired to change default ODE system settings
% dt = []; tspan = []; ode_params = {}; x0 = []; % change as desired to change default ODE system settings
ode_names = {'Linear','Logistic_Growth','Van_der_Pol','Duffing','Lotka_Volterra','Lorenz','Rossler','rational','Oregonator','Hindmarsh-Rose','Pendulum','custom'};
[weights,x,t,x0,ode_name,ode_params,rhs] = gen_data(ode_num,ode_params,tspan,x0,tol_ode);
[sigma,noise,xobs,tobs,noise_ratio_obs] = add_noise(x,noise_ratio,t);

%% Set parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% data smoothing
expected_coeff_ratio = 100;
[~,sigma_est] = findcornerpts(xobs(:,1),tobs);
smoothing_window = floor( sigma_est*expected_coeff_ratio );

%%% library
polys = 0:5;                        % monomial powers to include in library
trigs = [];                         % sine / cosine frequencies to include in library

%%% weak formulation
phi_class = 1;                      % 1 = piecewisepoly test functions, 2 = gaussian test functions
max_d = 1;                          % use {test function, derivative pair} {phi, -phi'} vs {phi', -phi''}, etc.
tau = 10^-16;                           % tau > 1: [tau,tauhat]=[m,p] directly. else...
tauhat = -2;                        % tauhat > 0: tauhat = width-at-half-max. tauhat < 0: use corner point
K_frac = 1;                         % K_frac > 1: K (number of test functions) = K_frac. K_frac < 1: K = T*K_frac (T = number of timespoints)
overlap_frac = 1;                   % fraction of allowable support overlap between neighboring test functions 
relax_AG =0;                        % convex combination between uniform and gradient-adapted placement of test functions

%%% optimization
scale_Theta = 2;                    % toggle normalize columns of Theta_0 by norm = scale_Theta_0
useGLS = 0;                         % useGLS = 0 calls OLS. useGLS > 0 calls GLS with cov (1-useGLS)*C + useGLS*diag(C), C = Vp*Vp'
lambda = 10.^(linspace(-4,0,100));  % lambda values for MSTLS
gamma = 10^-inf*sqrt(noise_ratio);                          % Tikhonov regularization
alpha_loss = 0.9;                   % convex combination between rel. resid. and rel. sparsity used in lambda loss function

%%% additional params for adaptive grid method (less important to tune)
overlap_frac_ag = 0.8;              % (involved in mini-routine to find steep gradients)
mt_ag_fac = [0;8]; pt_ag_fac=[0;2]; % (involved in mini-routine to find steep gradients)

%% Toggle Display results and/or Compare to SINDy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toggle_print_w = 2;                 % print weights
toggle_plot = 0;                    % print data, distribution of test functions, covariance sparsity
toggle_plot_ddd = 1;                % toggle simulate data-driven dynamics (ddd), 1 calls WSINDy weights, 2 calls SINDy weights, >2 displays error in vector field on a plot with toggle_plot_ddd points in each dimensions (dims 1,2,3) 
thresh = {'inf',3}; mult = 0;       % thresh - stop simulation when norm(x_dd,thresh{1})>thresh{2}*max(norm(x_true,thresh{1}) 
toggle_plot_resid = 0;              % plot residual with true weights and learned weights
toggle_plot_loss = 0;               % plot loss function
toggle_plot_derivs = 0;             % plot weak derivatives using adaptive grid (only for relax_AG>0)
toggle_plot_approx_sys = 0;         % plot trapezoidal-rule integral of Theta*weights 
toggle_plot_fft = 1;

%% get WSINDy model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

run_sindy = 1;                      % toggle get standard SINDy model
useFD = 0;                          % derivative method for standard SINDy.
                                    % useFD = centered finite difference stencil width (one side). useFD = 0 calls TVdiff
[w_sparse,w_sparse_sindy,true_nz_weights,loss_wsindy,loss_sindy,ET_wsindy,...
    ET_sindy,grids,pts,mts,Gs,bs,RTs,Theta_0,tags,bweaks,dxobs_0,vs] ...
    ...
    = wsindy_ode_fun(xobs,tobs,weights,...
    polys,trigs,...
    phi_class,max_d,tau,tauhat,K_frac,overlap_frac,relax_AG,...
    scale_Theta,useGLS,lambda,gamma,alpha_loss,...
    overlap_frac_ag,pt_ag_fac,mt_ag_fac,run_sindy,useFD,w);

%% Display results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rhs_dd = build_vector_field(w_sparse, tags);
[err_wsindy,err_sindy,tp_w,tp_s,t_dd,x_dd,F_dd,lambda_hat_w,lambda_hat_s,resid,resid_true] = display_results(w_sparse,true_nz_weights,w_sparse_sindy,loss_wsindy,...
    loss_sindy,lambda,noise_ratio,noise_ratio_obs,sigma,ET_wsindy,ET_sindy,xobs,x,tobs,t,grids,Gs,RTs,bs,Theta_0,...
    mts,pts,toggle_print_w,toggle_plot,toggle_plot_ddd,thresh,mult,toggle_plot_resid,...
    toggle_plot_loss,toggle_plot_derivs,toggle_plot_approx_sys,toggle_plot_fft,bweaks,useFD,dxobs_0,...
    tol_ode,x0,ode_name,ode_params,tags,vs);