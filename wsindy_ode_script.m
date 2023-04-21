%% Generate clean data (x,t,weights)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all; close all; clc;

ode_num = 6;                        % select ODE from list below
tol_ode = 1e-12;                    % set tolerance (abs and rel) of ode45

% tspan = [0:0.01:30]; ode_params = {0.1, 0.1, 5}; x0 = [0;1]; % example custom ODE system parameters
tspan = []; ode_params = {}; x0 = []; % ODE system parameters
ode_names = {'Linear','Logistic_Growth','Van_der_Pol','Duffing',... %1-4
             'Lotka_Volterra','Lorenz','Rossler','rational',...     %5-8
             'Oregonator','Hindmarsh-Rose','Pendulum','custom'};    %9-12
[weights,x,t,x0,ode_name,ode_params,rhs] = gen_data(ode_num,ode_params,tspan,x0,tol_ode);

%% load data from file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all;

% dr = '/home/danielmessenger/Dropbox/Boulder/research/data/WENDy_data/ode_data/'; 
% filename = 'gyroceptron_r_ep05.mat';
% load([dr,filename],'x','t','weights','ode_name','ode_params','features','params')

%% subsample data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subsamp = 1;                       % subsample data in time
tobs = t(1:subsamp:end); xsub = x(1:subsamp:end,:);
[M,nstates] = size(xsub);

%% add noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(1);
noise_ratio = 0.1;
noise_dist = 0;
noise_alg = 0;
rng_seed = rng().Seed; rng(rng_seed);
[xobs,noise,noise_ratio_obs,sigma] = gen_noise(xsub,noise_ratio,noise_dist,noise_alg);

%% Set parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% library
polys = [0:5];                        % Monomials. DEFAULT: polys = 0:5; 
trigs = [];                         % sine / cosine terms. DEFAULT: trigs = [];
filter = {};%{@(tags) tags(:,2)==0};
custom_tags = [];
custom_fcns = {};%prodlib(nstates,build_poly_lib(nstates,0:1,[],filter),build_trig_lib(nstates,[0 2 6],[],filter))';                    % custom terms. DEFAULT = {};
                                    % *(can be function handles or poly/trig tags) 

%%% weak formulation
phi_class = 1;%@(t)exp(9./(t.^2-1-eps)); % Choose test function class. DEFAULT: phi_class = 1. 
                                    % ---phi_class = 1: piecewisepoly test functions, 
                                    % ---phi_class = 2: gaussian test functions
                                    % ---phi_class = function_handle: use function handle and symbolically compute derivative
tau = 10^-16; tauhat = -1;          % Choose test function params. DEFAULT: [tau,tauhat] = [10^-16,-1]. 
                                    % ---tau > 1: [tau,tauhat] = [m,p] directly
                                    % ---tauhat > 0: tauhat = width-at-half-max. 
                                    % ---tauhat < 0: use corner point (recommended)
K_frac = 1000;                         % Choose number of test functions K. DEFAULT: K_frac = 1.
                                    % ---K_frac > 1: K = K_frac. 
                                    % ---K_frac <= 1: K = M*K_frac (M = number of timespoints)
%%% optimization
scale_Theta = 2;                    % rescale data. DEFAULT: scale_theta = 2.
                                    % ---scale_Theta < 0: rescale data: xobs -> xobs./(-scale_theta*rms(xobs))
                                    % ---scale_Theta > 0: rescale data: xobs -> xobs/||xobs||_{2p} where p=max(polys)
                                    % ---scale_theta == 0: no rescaling
lambda = 10.^(linspace(-4,0,100));  % sparsity factor(s). DEFAULT: lambda = 10.^(linspace(-4,0,100)).
alpha_loss = 0.8;                   % convex combination btw resid and sprasity in MSTLS loss. DEFAULT: 0.8.
                                    % *(smaller libraries require alpha_loss closer to 1)
gamma = 0;                          % Tikhonov regularization. DEFAULT: gamma = 0.
                                    % *(helps for ill-conditioned Theta, e.g. combined trig + poly libs)

%%% data smoothing 
smoothing_window = ceil(length(tobs)/100);  % rough guess for smoothing window. DEFAULT: smoothing_window = ceil(length(tobs)/100). 
                                    % *(should over-estimate the optimal window)
                                    % *(if smoothing not detected to be advantageous, window length will be set to 1)

%% additional params less important to tune
overlap_frac_ag = 0.8;              % (involved in mini-routine to find steep gradients)
mt_ag_fac = [1;0]; pt_ag_fac=[1;0]; % (involved in mini-routine to find steep gradients)
max_d = 1;                          % use {test function, derivative pair} {phi, -phi'} vs {phi', -phi''}, etc.
overlap_frac = 1;                   % fraction of allowable support overlap between neighboring test functions 
                                    % often the adaptive grid will pack
                                    % test functions very close together,
                                    % leading to nearly redundant rows.
                                    % restricting the support overlap prevents
                                    % this.
relax_AG = 0;                       % convex combination between uniform and gradient-adapted placement of test functions
                                    % adaptive grid helps for dynamics
                                    % with sharp transitions (e.g. Van der
                                    % Pol)
useGLS = 0;                     % useGLS = 0 calls OLS. useGLS > 0 calls GLS with cov (1-useGLS)*C + useGLS*diag(C), C = Vp*Vp'
                                    % GLS is helpful when the test
                                    % functions have very small support and
                                    % when the jacobian of the true system
                                    % is small (e.g. linear Dynamics)

%% Toggle Display results and/or Compare to SINDy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toggle_print_w = 2;                 % print weights
toggle_plot = 0;                    % print data, distribution of test functions, covariance sparsity
toggle_plot_ddd = 1;                % toggle simulate data-driven dynamics (ddd), 1 calls WSINDy weights, 2 calls SINDy weights, >2 displays error in vector field on a plot with toggle_plot_ddd points in each dimensions (dims 1,2,3) 
thresh = {'inf',3}; mult = 0;       % thresh - stop simulation when norm(x_dd,thresh{1})>thresh{2}*max(norm(x_true,thresh{1}) 
toggle_plot_resid = 1;              % plot residual with true weights and learned weights
toggle_plot_loss = 1;               % plot loss function
toggle_plot_derivs = 0;             % plot weak derivatives using adaptive grid (only for relax_AG>0)
toggle_plot_approx_sys = 0;         % plot trapezoidal-rule integral of Theta*weights 
toggle_plot_fft = 1;                % plot data power spectrum
toggle_plot_filter_weights = 1;     % plot smoothing filters

%% toggle run SINDy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

run_sindy = 0;                      % toggle get standard SINDy model
useFD = 2;                          % derivative method for standard SINDy. DEFAULT: useFD = 0.
                                    % ---useFD > 0: use centered finite difference w/ width 2*useFD+1
                                    % ---useFD == 0: use TVdiff with reg_param = 1/sqrt(M), M=num timepoints
                                    % ---useFD < 0: use TVdiff with reg_param = (-useFD)/sqrt(M)

 %% get WSINDy model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                   
[w_sparse,w_sparse_sindy,true_nz_weights,loss_wsindy,loss_sindy,ET_wsindy,...
    ET_sindy,grids,pts,mts,Gs,bs,M_diag,RTs,Theta_0,tags,bweaks,dxobs_0,vs,filter_weights,xobs_smooth] ...
    = wsindy_ode_fun(xobs,tobs,weights,...
    polys,trigs,custom_tags,custom_fcns,...
    phi_class,max_d,tau,tauhat,K_frac,overlap_frac,relax_AG,...
    scale_Theta,useGLS,lambda,gamma,alpha_loss,...
    overlap_frac_ag,pt_ag_fac,mt_ag_fac,run_sindy,useFD,smoothing_window);

%% Display results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tol_ode_FS = 10^-10;
x0_FS = x(1,:);

[err_wsindy,err_sindy,tp_w,tp_s,t_dd,x_dd,F_dd,lambda_hat_w,lambda_hat_s,resid,resid_true,x_approx] = ...
    display_results(w_sparse,true_nz_weights,w_sparse_sindy,loss_wsindy,...
    loss_sindy,lambda,noise_ratio,noise_ratio_obs,sigma,ET_wsindy,ET_sindy,xobs_smooth,x,tobs,t,grids,Gs,RTs,bs,Theta_0,M_diag,...
    mts,pts,toggle_print_w,toggle_plot,toggle_plot_ddd,thresh,mult,toggle_plot_resid,...
    toggle_plot_loss,toggle_plot_derivs,toggle_plot_approx_sys,toggle_plot_fft,toggle_plot_filter_weights,bweaks,useFD,dxobs_0,...
    tol_ode_FS,x0_FS,ode_name,ode_params,tags,vs,filter_weights);