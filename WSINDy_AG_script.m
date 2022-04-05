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

dt = 0.01;

ode_num = 7;                        % select ODE system from the list ode_names
tol_ode = 1e-12;                    % ode45 tolerance (abs and rel) for generating data
num_runs = 1;                       % run script this many times

noise_ratio = 0;                    % set signal-to-noise ratio (L2 sense)
rng('shuffle');                   % comment out to use same noise as previous run for reproducibility
rng_seed =rng; 
rng_seed =rng_seed.Seed;
rng(rng_seed);

useFD_SINDy = 2;                          % SINDy finite difference differentiation order, if =0, then uses TVdiff

use_preset_params = 0;              % use parameters in script gen_data.m
if ~use_preset_params
    % library and optimization params
    polys = 0:4;                    % polynomial degrees to include
    trigs = [];                     % trig frequencies to include
    gamma = 10^-Inf;                % Tikhonov regularization
    lambda_mult = 4;                % sparsity param = min_{w*_i~=0}|w*_i|/lambda_mult
    scale_Theta = 0;                % normalize columns of Theta (needs work)

    % set WSINDy params
    tau_p = 15;                      % sets poly degree p = -tau_p. If tau_p<0, instead sets p by enforcing that test function has value 10^-tau_p at penultimate support point 
    r_whm = 0.1/dt;                     % test function whm - sets support length
    tau = 0;                        % toggle between uniform and adapted grid
    K = 1000;                        % number of test functions per coordinate (exact for uniform, approximate for adaptive)
    p = 2; s = 16;                  % (parameters for adaptive grid)
    useGLS = 0;                     % toggle use generalized least squares
end


toggle_display_weights = 2;
toggle_plot = 1;

% ode params
ode_names = {'Linear','Logistic_Growth','Van_der_Pol','Duffing','Lotka_Volterra','Lorenz','Custom'};
ode_name = ode_names{ode_num};
if strcmp(ode_name, 'Linear') %1
    ode_params = {[[-0.1 2];[-2 -0.1]]};x0 = [3;0]; tspan = 0:dt:15;

elseif strcmp(ode_name, 'Logistic_Growth') %2
    ode_params = {2}; x0 = 0.01; tspan = 0:0.005:10;

elseif strcmp(ode_name, 'Van_der_Pol') %3
    ode_params = {4}; x0 = [0;1]; tspan = 0:dt:30;

elseif strcmp(ode_name, 'Duffing') %4
    mu = 0.2;ode_params = {mu, mu^2/4*5,1}; x0 = [0;2]; tspan = 0:0.01:30;       

elseif strcmp(ode_name, 'Lotka_Volterra') %5
    alpha=2/3;beta = 4/3; ode_params = {alpha,beta,1,1}; x0 = [10;10]; tspan = 0:0.02:200;

elseif strcmp(ode_name, 'Lorenz') %6
    ode_params = {10, 8/3, 28}; tspan = .001:dt:5; x0 = [-8 7 27]';
    %x0 = [rand(2,1)*30-15;rand*30+10]
elseif strcmp(ode_name, 'Custom') %7
    mu=-1;
    tspan = .01:dt:30; x0 = [0.3 0.3]';
    ode_func = @(x) [x(1)-x(1).^2;-0.1*x(2)+0.5*x(1)+0.5];
    weights = {[[1 0 1];[2 0 -1]], [[0 1 -0.1];[1 0 0.5];[0 0 0.5]]};
    ode_params = {ode_func, weights}; 
    %x0 = [rand(2,1)*30-15;rand*30+10]
end

tps_list = zeros(num_runs,2);
errs_list = zeros(num_runs,4);

for nn=1:num_runs

    gen_data;

    %%% get WSINDy model

    get_wsindy_model;
    tpr_wsindy = tpscore(w_sparse,true_nz_weights);
    err_wsindy = [norm(true_nz_weights(:)-w_sparse(:))/norm(true_nz_weights(:)) ...
        norm(true_nz_weights(true_nz_weights~=0)-w_sparse(true_nz_weights~=0))/norm(true_nz_weights(:))];
    tps_list(nn,1) = tpr_wsindy;
    errs_list(nn,1:2) = err_wsindy;
    
    %%% get SINDy model

    if useFD_SINDy>=0
        tic;
        w_sparse_sindy = standard_sindy(t,xobs,Theta_0,M_diag, useFD_SINDy,n,lambda,gamma);
        ET_s = toc;
        err_sindy = [norm(w_sparse_sindy(:)-true_nz_weights(:));norm(w_sparse_sindy(true_nz_weights~=0)-true_nz_weights(true_nz_weights~=0))]/norm(true_nz_weights(:));
        tpr_sindy = tpscore(w_sparse_sindy,true_nz_weights);
        tps_list(nn,2) = tpr_sindy;
        errs_list(nn,3:4) = err_sindy;
    else
        w_sparse_sindy = w_sparse*0;
        ET_s = 0;
        err_sindy = [NaN NaN];
        tpr_sindy = NaN;
        tps_list(nn,2) = tpr_sindy;
        errs_list(nn,3:4) = err_sindy;
    end        
    
    %%% visualize basis and covariance, display error analysis in command window

%     figure(3); clf    set(gcf, 'position', [1250 10 700 450])
    
    if toggle_plot
        num_rows = 2+double(useGLS>0);
        for d=1:n
            subplot(num_rows,n,d) 
            plot(tobs,xobs(:,d),'r-',tobs(floor(mean(ts_grids{d},2))),mean(xobs(:,d))*ones(size(ts_grids{d})),'.k')
            subplot(num_rows,n,n+d)
            plot(tobs,mats{d}{1}')
            if num_rows==3
                subplot(num_rows,n,2*n+d)
                spy(RTs{d})
            end
        end
        drawnow
    end
    
    supp_range = [];
    for i=1:n    
        supp_range = [supp_range; ts_grids{i}*[-1;1]];
    end

    clc;
    disp(['log10 2norm err (all weights) (WSINDy)=',num2str(log10(err_wsindy(1)))])
    disp(['log10 2norm err (all weights) (SINDy)=',num2str(log10(err_sindy(1)))])
    disp(['log10 2norm err (true nz weights) (WSINDy)=',num2str(log10(err_wsindy(2)))])
    disp(['log10 2norm err (true nz weights) (SINDy)=',num2str(log10(err_sindy(2)))])
    disp(['----------------'])
    disp(['TPR (WSINDy)=',num2str(tpr_wsindy)])
    disp(['num times identified (WSINDy) =',num2str(length(find(tps_list(:,1)==1))),'/',num2str(nn)])
    disp(['avg err (WSINDy) =',num2str(mean(errs_list(1:nn,1)))])
    disp(['TPR (SINDy)=',num2str(tpr_sindy)])
    disp(['num times identified (SINDy) =',num2str(length(find(tps_list(:,2)==1))),'/',num2str(nn)])
    disp(['avg err (SINDy) =',num2str(mean(errs_list(1:nn,3)))])
    disp(['----------------'])
    disp(['min/mean/max deg=',num2str([min(ps_all) mean(ps_all) max(ps_all)])])
    disp(['min/mean/max supp=',num2str([min(supp_range) mean(supp_range) max(supp_range)])])
    disp(['Num test fcns=',num2str(size(V,1))])
    disp(['Num trial fcns=',num2str(size(G,2))])
    disp(['----------------'])
    disp(['Run time (WSINDy) =',num2str(ET)])
    disp(['Run time (SINDy) =',num2str(ET_s)])
    
    if toggle_display_weights == 1
        disp(['WSINDy, SINDy, true, library:'])
        disp([w_sparse w_sparse_sindy true_nz_weights tags])
    elseif toggle_display_weights == 2
        disp(['WSINDy, true, lib:'])
        disp([w_sparse true_nz_weights tags])
    elseif toggle_display_weights == 3
        disp(['WSINDy, true:'])
        disp([w_sparse true_nz_weights])
    end
    
end