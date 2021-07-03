%% Generate Clean Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ode_num = 3;                         % select ODE from list below. 1,3,4,5,6,11 used in manuscript 
tol_ode = 1e-10;                     % set tolerance (abs and rel) of ode45

ode_names = {'Linear','Logistic_Growth','Van_der_Pol','Duffing',...
    'Lotka_Volterra','Lorenz','Rossler','rational','Oregonator',...
    'Hindmarsh-Rose','Pendulum','custom'};

tspan = []; ode_params = []; x0 = []; % manually input timepoints, ODE parameters, initial conditions   
[weights,x,t,x0,ode_name,ode_params] = gen_data(ode_num,ode_params,tspan,x0,tol_ode);

%% Set WSINDy Hyperparameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% library hyperparameters
polys = 0:8;                        % monomial powers to include in library
trigs = [];                      % trig frequencies to include in library (includes sin and cos separately)

% test function hyperparameters: phi(t) = (1-(t/(m_t*dt))^2)^(p_t)
mt = 50;                            % 2*m_t+1 = number of timepoints that test function is supported on
pt = 5;                             % test function polynomial degree
st = 10;                   % gap between test function centers

% optimization hyperparameters
lambda = 0.25;                       % sparsity threshold
gamma = 0.00;                          % Tikhonov regularization, so that weights w are chosen by minimizing ||G*w-b||^2+gamma^2*||w||^2
    
%% Add Noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

noise_ratio = 0.1;                 % set ratio of noise to signal (L2 sense), if negative, sets sigma to -noise_ratio
rng('shuffle');                     % comment out to reproduce previous noise 
rng_seed =rng().Seed;
rng(rng_seed);
[sigma,noise,xobs,tobs,noise_ratio_obs] = add_noise(x,noise_ratio,t);

%% Build Dictionary and run WSINDy 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic,
% generate true coefficients if available
[m,n] = size(xobs); 
tags = get_tags(polys,trigs,n);
J = size(tags,1);
if ~isempty(weights)
    true_nz_weights = get_true_weights(weights,tags,n);
else
    true_nz_weights = ones(J,n)/J;
end

% build Theta matrix    
Theta_0 = build_theta(xobs,tags);

% get test function values on timegrid
dt = mean(diff(tobs));
Cfs = phi_int_weights(mt,pt);
v = Cfs(end-1,:)*dt;
vp = Cfs(end,:)*(mt*dt)^(-1)*dt;
grid_tf = 1:st:m-2*mt;

% build linear system
G = conv2(v,1,Theta_0,'valid');
G = G(grid_tf,:);
b = conv2(vp,1,xobs,'valid');
b = b(grid_tf,:);

% get model ceofficients
[w_wsindy,its_w] = sparsifyDynamics(G,b,lambda,gamma);
ET_wsindy = toc;

%% Run standard SINDy
run_standard_sindy = 1;             % =1, run standard SINDy; note all timepoints are used 
useFD = 1;                          % Centered finite difference stencil width (one side). =0 calls TVdiff 

if run_standard_sindy
tic;
    if useFD<=0
        L = length(tobs)-2;
        C = fdcoeffF(1, tobs(2), tobs(1:3));
        D = spdiags(repmat(C(:,end)',L,1), 0:2 ,L,length(tobs));
        naive_diff = D*xobs;
        dxobs_0 = 0*naive_diff;
        if useFD ==0
            reg_param = dt;
        else
            reg_param = (-useFD)*dt;
        end
        for j=1:n
            dxobs_0(:,j) = TVRegDiff( xobs(2:end-2,j), 20, reg_param, naive_diff(:,j), [], [], dt, 0, 0 );
        end
        Theta_sindy = Theta_0(max(useFD,1)+1:end-max(useFD,1),:);
    else
        L = length(tobs)-useFD*2;
        C = fdcoeffF(1, tobs(useFD+1), tobs(1:useFD*2+1));
        D = spdiags(repmat(C(:,end)',L,1), 0:useFD*2 ,L,length(tobs));
        dxobs_0 = D*xobs;
        Theta_sindy = Theta_0(max(useFD,1)+1:end-max(useFD,1),:);    
    end
[w_sindy,its_s] = sparsifyDynamics(Theta_sindy,dxobs_0,lambda,gamma);
ET_sindy = toc;
else
    w_sindy = 0*w_wsindy;
    ET_sindy = 0;
    its_s = 0;
end
    
%% Toggle Display results and/or Compare to SINDy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toggle_print_w = 1;                 % print weights
toggle_plot = 1;                    % print data, distribution of test functions, covariance sparsity
toggle_plot_ddd = 1;                % toggle simulate data-driven dynamics (ddd), 1 calls WSINDy weights, 2 calls SINDy weights, >2 displays error in vector field on a plot with toggle_plot_ddd points in each dimensions (dims 1,2,3) 
thresh = {'inf',3}; mult = 0;       % thresh - stop simulation when norm(x_dd,thresh{1})>thresh{2}*max(norm(x_true,thresh{1}) 
toggle_plot_resid = 1;              % plot residual with true weights and learned weights

clc;
close all; 

err_wsindy = [norm(w_wsindy(:)-true_nz_weights(:));norm(w_wsindy(true_nz_weights~=0)-true_nz_weights(true_nz_weights~=0))]/norm(true_nz_weights(:));
err_sindy = [norm(w_sindy(:)-true_nz_weights(:));norm(w_sindy(true_nz_weights~=0)-true_nz_weights(true_nz_weights~=0))]/norm(true_nz_weights(:));        

disp(['log10 2norm err (all weights) (WSINDy)=',num2str(log10(err_wsindy(1)))])
disp(['log10 2norm err (all weights) (SINDy)=',num2str(log10(err_sindy(1)))])
disp(['log10 2norm err (true nz weights) (WSINDy)=',num2str(log10(err_wsindy(2)))])
disp(['log10 2norm err (true nz weights) (SINDy)=',num2str(log10(err_sindy(2)))])
tp_w = tpscore(w_wsindy,true_nz_weights);
tp_s = tpscore(w_sindy,true_nz_weights);
disp(['TPR (WSINDy, SINDy)=',num2str([tp_w tp_s])])
disp(' ')
disp(['Noise_ratio, sigma =',num2str([noise_ratio_obs sigma])])
disp(['Run time ([WSINDy SINDy]) =', num2str([ET_wsindy ET_sindy])])
disp(['STLS its (WSINDy, SINDy)=', num2str([its_w its_s])])
disp(' ')
disp(['Num timepoints =',num2str(size(xobs,1))])
disp(['Num Trial Fcns =',num2str(size(w_wsindy,1))])
disp(['Num Test Fcns  =',num2str(length(grid_tf))])
disp(['Test function degree  =',num2str(pt)])
disp(['Test function support =',num2str(2*mt+1)])
disp(['log10(cond([G Theta_0])) =',num2str(cellfun(@(x) log10(cond(x)),[{G} {Theta_0}]))])
disp(' ')

if toggle_print_w == 1
    disp(['true, wsindy, sindy'])
    disp([true_nz_weights w_wsindy w_sindy])
end

if toggle_plot>0
    figure(1); clf
    set(gcf, 'units','normalized','outerposition',[0 0.5 0.5 0.5])
    for nn=1:n
        subplot(1,n+1,nn) 
        plot(tobs,xobs(:,nn),'r-',tobs(grid_tf+mt),mean(xobs(:,nn))*ones(length(grid_tf),1),'.k')
        legend({['xobs_{',num2str(nn),'}'],'test function centers'})
    end
    subplot(1,n+1,n+1) 
    plot(tobs(1:2*mt+1),Cfs(1,:),'b-',tobs(1:2*mt+1),Cfs(2,:),'g-')
    legend({['\phi'],'\phi^\prime'})
end

if or(toggle_plot_ddd==1,toggle_plot_ddd==2)
    figure(2); clf
    set(gcf, 'units','normalized','outerposition',[0.5 0.5 0.5 0.5])
    if toggle_plot_ddd == 1
        w = w_wsindy;
    elseif toggle_plot_ddd == 2
        w = w_sindy;
    end
    [~,~,t_dd,x_dd,F_dd] = view_ddd_fcn(thresh,mult,w,tol_ode,x0,tobs(1:5:end),tol_ode,ode_name,ode_params,tags,toggle_plot_ddd);
    title(num2str(tobs(end)))
else
    t_dd = [];
    x_dd = [];
    F_dd = [];
end

if toggle_plot_resid
    figure(3); clf
    set(gcf, 'units','normalized','outerposition',[0.5 0 0.5 0.5])
    grids_r = grid_tf+mt;
    resid = (b-G*w_wsindy)./vecnorm(b);
    resid_true = (b-G*true_nz_weights)./vecnorm(b);
    ylims = max(noise_ratio,max(abs(resid(:)))); ylims = [-ylims ylims];
    for nn=1:n
        subplot(2,n,nn)
        plot(t,x(:,nn)*(ylims(2)/max(x(:,nn))),'r--',t(grids_r),resid(:,nn),'k.')
        ylim(ylims)
        legend({'clean data (for reference)','WSINDy resid'})
        ylabel(['x_{',num2str(nn),'}'])
        subplot(2,n,n+nn)        
        plot(t,x(:,nn)*(ylims(2)/max(x(:,nn))),'r--',t(grids_r),resid_true(:,nn),'k.')
        ylim(ylims)
        legend({'clean data (for reference)','True resid'})
        ylabel(['x_{',num2str(nn),'}'])
    end
else
    resid = [];
    resid_true = [];
end