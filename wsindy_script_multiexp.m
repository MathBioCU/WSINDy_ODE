%%% choose ode
ode_names = {'Linear','Logistic_Growth','Van_der_Pol','Duffing',... %1-4
             'Lotka_Volterra','Lorenz','Rossler','rational',...     %5-8
             'Oregonator','Hindmarsh-Rose','Pendulum','custom', 'G_bhatia_1','G_bhatia_2','G_bhatia_3'};    %9-15
ode_num = 14;                        % select ODE from list below
ode_params = {}; % fix ode params across experiments

% number of experiments
num_exp = 10; 
x_cell = cell(num_exp,1); 

% choose initial conditions for each exp
x0s = linspace(0.1,0.7,num_exp)';
x0s = [x0s flipud(x0s)];
% x0s = rand(num_exp,2);
nstates = size(x0s,2);

% choose timespan for each exp
t_cell = repmat({[0:0.05:20]},num_exp,1); %%% for each initial condition, generate data

% gen data
tol_ode = 1e-12;                    % set tolerance (abs and rel) of ode45
for j=1:num_exp
    [weights,x_cell{j},t_cell{j},~,ode_name,ode_params,~] = gen_data(ode_num,ode_params,t_cell{j},x0s(j,:),tol_ode);
end

%%% load data from file
% load([dr,filename],'x_cell','tspans')


%%% subsample data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subsamp = 1;                       % subsample data in time
tobs_cell = cellfun(@(t)t(1:subsamp:end),t_cell,'uni',0);
xsub_cell = cellfun(@(x)x(1:subsamp:end,:),x_cell,'uni',0);

%% add noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng('shuffle');
noise_ratio = 0.0;
noise_dist = 0;
noise_alg = 0;
rng_seed = rng().Seed; rng(rng_seed);
xobs_cell = xsub_cell;
for j=1:num_exp 
    [xobs_cell{j},noise,noise_ratio_obs,sigma] = gen_noise(xsub_cell{j},noise_ratio,noise_dist,noise_alg);
end

%% Set parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% library
polys = [0:2];                        % Monomials. DEFAULT: polys = 0:5; 
trigs = [];                         % sine / cosine terms. DEFAULT: trigs = [];
filter = {};%{@(tags) tags(:,2)==0};
custom_tags = [];
custom_fcns = {};%prodlib(nstates,build_poly_lib(nstates,0:1,[],filter),build_trig_lib(nstates,[0 2 6],[],filter))';                    % custom terms. DEFAULT = {};
                                    % *(can be function handles or poly/trig tags) 

%%% weak formulation
phi_class = @(t)exp(9./(t.^2-1-eps)); % Choose test function class. DEFAULT: phi_class = 1. 
                                    % ---phi_class = 1: piecewisepoly test functions, 
                                    % ---phi_class = 2: gaussian test functions
                                    % ---phi_class = function_handle: use function handle and symbolically compute derivative
tau = 10^-16; tauhat = -2;          % Choose test function params. DEFAULT: [tau,tauhat] = [10^-16,-1]. 
% tau = 15; tauhat = 5;          % Choose test function params. DEFAULT: [tau,tauhat] = [10^-16,-1]. 

                                    % ---tau > 1: [tau,tauhat] = [m,p] directly
                                    % ---tauhat > 0: tauhat = width-at-half-max. 
                                    % ---tauhat < 0: use corner point (recommended)
K_frac = 1000;                         % Choose number of test functions K. DEFAULT: K_frac = 1.
                                    % ---K_frac > 1: K = K_frac. 
                                    % ---K_frac <= 1: K = M*K_frac (M = number of timespoints)
%%% optimization
scale_Theta = 0;                    % rescale data. DEFAULT: scale_theta = 2.
                                    % ---scale_Theta < 0: rescale data: xobs -> xobs./(-scale_theta*rms(xobs))
                                    % ---scale_Theta > 0: rescale data: xobs -> xobs/||xobs||_{2p} where p=max(polys)
                                    % ---scale_theta == 0: no rescaling
lambda = 10.^(linspace(-5,0,100));  % sparsity factor(s). DEFAULT: lambda = 10.^(linspace(-4,0,100)).
alph = 0.01;                   % convex combination btw resid and sprasity in MSTLS loss. DEFAULT: 0.8.
                                   % *(smaller libraries require alpha_loss closer to 1)
gamma = 0;                          % Tikhonov regularization. DEFAULT: gamma = 0.
                                    % *(helps for ill-conditioned Theta, e.g. combined trig + poly libs)

%%% data smoothing 
smoothing_window = 0;  % rough guess for smoothing window. DEFAULT: smoothing_window = ceil(length(tobs)/100). 
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

%% WSINDy multiexp
tic,
%%% get linear system
                             
xobs_filter = xobs_cell;
G_cell = cell(num_exp,nstates);
b_cell = cell(num_exp,nstates);
mt_cell = cell(num_exp,1);
for j=1:num_exp
[grids,pts,mt_cell{j},G_cell(j,:),b_cell(j,:),M_diag,RTs,Theta,tags,bweaks,vs,filter_weights,xobs_filter{j}] ...
    = get_Gb(xobs_cell{j},tobs_cell{j},polys,trigs,custom_tags,custom_fcns,...
        phi_class,max_d,tau,tauhat,K_frac,overlap_frac,relax_AG,...
        scale_Theta,useGLS,overlap_frac_ag,pt_ag_fac,mt_ag_fac,smoothing_window);
end
Gs = arrayfun(@(i)cell2mat(G_cell(:,i)),1:nstates,'uni',0);
bs = arrayfun(@(i)cell2mat(b_cell(:,i)),1:nstates,'uni',0);

%%% regress
loss_wsindy = cell(nstates,1);
w_sparse = cell(1,nstates);
for nn=1:nstates
    M_diag = ones(size(Gs{nn},2),1);
    if length(lambda)==1
        [w_sparse{nn},~] = sparsifyDynamics(Gs{nn},bs{nn},lambda,1,gamma,M_diag);
        loss_wsindy = [];
    else
        if ~isempty(M_diag)
            M_scale_b = [1;M_diag];
        else
            M_scale_b = [];
        end
        alpha_loss = 1/(alph*size(Gs{nn},2)+1);
        [w_sparse{nn},loss_wsindy{nn},~] = wsindy_pde_RGLS_seq(lambda,gamma,[bs{nn} Gs{nn}],1,M_scale_b,alpha_loss);
    end
end
w_sparse = cell2mat(w_sparse);
ET_wsindy = toc;

%% Toggle Display results and/or Compare to SINDy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Display results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

exp_test = num_exp;
toggle_print_w = 2;                 % print weights
toggle_plot_ddd = 1;                % toggle simulate data-driven dynamics (ddd), 1 calls WSINDy weights, 2 calls SINDy weights, >2 displays error in vector field on a plot with toggle_plot_ddd points in each dimensions (dims 1,2,3) 
thresh = {'inf',3}; mult = 0;       % thresh - stop simulation when norm(x_dd,thresh{1})>thresh{2}*max(norm(x_true,thresh{1}) 
toggle_plot_loss = 1;               % plot loss function
toggle_plot_fft = 1;                % plot data power spectrum
toggle_plot_filter_weights = 1;     % plot smoothing filters

true_nz_weights = get_true_weights(weights,tags,nstates);

G_0 = blkdiag(Gs{:});
b_0 = cell2mat(bs(:));

res = norm(G_0*reshape(w_sparse./M_diag,[],1)-b_0)/norm(b_0);
disp(['rel. resid (WSINDy)=',num2str(res)])
if ~isempty(true_nz_weights)
    err_wsindy = [norm(w_sparse(:)-true_nz_weights(:));norm(w_sparse(true_nz_weights~=0)-true_nz_weights(true_nz_weights~=0))]/norm(true_nz_weights(:));
    disp(['log10 2norm err (all weights) (WSINDy)=',num2str(log10(err_wsindy(1)))])
    disp(['log10 2norm err (true nz weights) (WSINDy)=',num2str(log10(err_wsindy(2)))])
    tp_w = tpscore(w_sparse,true_nz_weights);
    disp(['TPR =',num2str([tp_w])])
    disp(' ')
else
    err_wsindy=NaN;err_sindy = NaN;
    tp_w=NaN; tp_s=NaN;
end

lambda_loss_wsindy = [];
for nn=1:nstates
    if length(lambda)>1
        indtemp = find(loss_wsindy{nn}(4,:)==0,1)-1; if isempty(indtemp);indtemp=length(lambda);end
        lambda_loss_wsindy(nn,:) = [indtemp lambda(indtemp)];
    else
        lambda_loss_wsindy(nn,:) = [1 lambda];
    end
end

disp(['Noise_ratio, sigma =',num2str([noise_ratio_obs sigma])])
disp(['Run time (WSINDy) =',num2str(ET_wsindy)])
lambda_hat_w  = lambda_loss_wsindy(:,2)';
disp(['lambda_hat_w  =', num2str(lambda_hat_w)])
disp(' ')
disp(['Num Trial Fcns =',num2str(size(w_sparse,1))])
disp(['Num Test Fcns =',num2str(cellfun(@(x) size(x,1),Gs))])
disp(['log10(cond(G)) =',num2str(cellfun(@(x) log10(cond(x)),Gs))])
disp(' ')

if and(toggle_print_w, ~isempty(true_nz_weights))
    disp(['tags true, wsindy'])
    disp([tags zeros(size(tags,1),1) true_nz_weights zeros(size(tags,1),1) w_sparse])
end

tol_ode_FS = 10^-10;
x0_FS = xsub_cell{exp_test}(1,:);
[M,nstates] = size(xsub_cell{exp_test});

xobs = xobs_cell{exp_test};
tobs = tobs_cell{exp_test};

if toggle_plot_ddd
    figure(2); clf
    w = w_sparse;
    x0_new=x0_FS+mult*(rand(size(x0_FS))-0.5);    
    [~,t_test,x_test,~] = sim_ode(x0_new,tobs,tol_ode,ode_name,ode_params,'s');
    [t_dd,x_dd,F_dd] = generate_ddd(w,tags,tol_ode,x0_new,tobs,norm(x_test,thresh{1})*thresh{2});
    
    if nstates>1
        for d=1:nstates
            legs = {};
            subplot(2*nstates,1,d) 
            plot(t_test,x_test(:,d),'k','linewidth',2);  legs{1} = 'true';
            hold on 
            plot(t_dd,x_dd(:,d),'b--','linewidth',2); legs{2} = 'learned';
            hold off
            try
                title(['rel L2 err=',num2str(norm(x_dd(:,d)-x_test(:,d))/norm(x_test(:,d)))])
            catch
                disp(['data-driven simulation halted early'])
            end
            legend(legs)
        end    
        if nstates==2
            legs = {};
            subplot(2*nstates,1,nstates+1:2*nstates)
            plot(x_test(:,1),x_test(:,2),'k','linewidth',2);  legs{1} = 'true';
            hold on 
            plot(x_dd(:,1),x_dd(:,2),'b--','linewidth',2); legs{2} = 'learned';
            hold off
        elseif nstates==3
            legs = {};
            subplot(2*nstates,1,nstates+1:2*nstates)
            plot3(x_test(:,1),x_test(:,2),x_test(:,3),'k','linewidth',2);  legs{1} = 'true';
            hold on 
            plot3(x_dd(:,1),x_dd(:,2),x_dd(:,3),'g--','linewidth',2); legs{2} = 'learned';
            hold off
        end
    else
        plot(t_dd,x_dd,'g-'); legs{1} = 'DD';
        hold on 
        plot(t_test,x_test,'k');  legs{2} = 'true';
        hold off
        legend(legs)
    end
    sgtitle('data-driven dynamics')
end
    
if and(length(lambda)>1,toggle_plot_loss)
    figure(4); clf
    for nn=1:nstates
        subplot(nstates,1,nn)
        semilogx(lambda,loss_wsindy{nn}(1,:),'ob-',lambda_loss_wsindy(nn,2),loss_wsindy{nn}(3,lambda_loss_wsindy(nn,1)),'rx'); hold on;
        legend({'wsindy',''}, 'location', 'best')
    end
    sgtitle('MSTLS loss function')
end

if toggle_plot_fft>0
    figure(7);clf
    for j=1:nstates
        m = (length(vs{j,1})-1)/2;
        Ufft = abs(fft(xobs(:,j)));
        Ufft = Ufft(floor(end/2):end);
        L = length(Ufft)-1;
        ks = -L:L;
        Ufft = [Ufft; flipud(Ufft(1:end-1))]/max(Ufft);
        subplot(nstates,1,j)
            semilogy(ks,Ufft)
            hold on
            Cfs_ffts = fft([zeros(1,M-2*m-1) vs{j,1}]);
            Cfs_ffts = abs(Cfs_ffts(floor(end/2):end));
            Cfs_ffts = [Cfs_ffts fliplr(Cfs_ffts(1:end-1))];
            Cfs_ffts = Cfs_ffts/max(Cfs_ffts);
            semilogy(ks,Cfs_ffts)
            [corner,~] = findcornerpts(xobs(:,j),tobs);
            k = corner(2);
            semilogy([-k k],Ufft(max(L+1-k,1))*[1 1],'o','markersize',12)
            hold off     
            ylim([min(Ufft)*0.1 max(Ufft)])
            xlim([min(ks) max(ks)])
            legend({'$\mathcal{F}(y)$','$\mathcal{F}(\phi)$','$k^*$'},'interpreter','latex','fontsize',14)
            title(['coord ',num2str(j)])
    end
    xlabel('k')
    sgtitle('FFT of data and test functions')
end

if and(toggle_plot_filter_weights>0,~isempty(filter_weights))
    legs = {};
    figure(8); clf
    for j=1:nstates
        m = (length(filter_weights{j})-1)/2;
        x = -m:m;
        plot(x,filter_weights{j},'o-','linewidth',2);
        hold on
        legs{j} = (['coord ',num2str(j)]);
    end
    ylim([0 max(cellfun(@(x)max(x),filter_weights))])
    legend(legs)
    title('filter weights')
    hold off;
end

function [grids,pts,mts,Gs,bs,M_diag,RTs,Theta,tags,bweaks,vs,filter_weights,xobs] ...
    = get_Gb(xobs,tobs,polys,trigs,custom_tags,custom_fcns,...
    phi_class,max_d,tau,tauhat,K_frac,overlap_frac,relax_AG,...
    scale_Theta,useGLS,overlap_frac_ag,pt_ag_fac,mt_ag_fac,smoothing_window)

    %%% get true weight vector
    n = size(xobs,2);
    m = length(tobs);
    tags = get_tags(polys,trigs,n);
    tags = unique([tags;custom_tags],'rows');
    J = size(tags,1);

    %%% get scales
    if scale_Theta == 0
        scale_x = [];
    elseif scale_Theta < 0
        scale_x = rms(xobs)*(-scale_Theta);
    else
        if ~isempty(real(tags))
            scale_x = (vecnorm(xobs.^max(real(tags)),2)).^(1./max(real(tags)));
        else
            scale_x =[];
        end
    end
   
    %%% set weak form quantities
    tic,

    %%% set number of test functions. *** decreases when overlapfrac<1
    if K_frac>1
        K = ceil(K_frac);
    elseif K_frac>0
        K = floor(m*K_frac);
    elseif K_frac<0
        K = J*ceil(-K_frac);
    end

    %%% set mts, pts
    if any(tau>1)
        if length(tau)==n
            mts = tau;
        else
            mts = tau*ones(n,1);
        end
        if length(tauhat)==n
            pts = tauhat;
        else
            pts = tauhat*ones(n,1);
        end
    else
        if tauhat > 0
            [pts,mts] = test_fcn_param_whm(tauhat,tobs,tau);
            pts = repmat(pts,n,1);
            mts = repmat(mts,n,1);
        elseif tauhat < 0  
            [mts,pts,~,~] = findcorners(xobs,tobs,tau,-tauhat,phi_class);
        else
            disp('Error: tauhat cannot be 0')
            return
        end
    end

        %%% smooth data
    if smoothing_window>0
%         smoothing_grid = linspace(-1,1,2*w+3);
%         filter_weights = 1-smoothing_grid.^2;
%         filter_weights = filter_weights/sum(filter_weights);
%         filter_weights = filter_weights(2:end-1);
%         xobs = [flipud(xobs(2:w+1,:));xobs;flipud(xobs(end-w:end-1,:))];
%         xobs = conv2(filter_weights,1,xobs,'valid');
        filter_weights = cell(n,1);
        for j=1:n
%             filter_weights{j} = get_smoothing_weights(xobs(:,j),tobs,smoothing_window);
            [filter_weights{j},~,~] = get_optimal_SMAF(tobs,xobs(:,j),[],smoothing_window,[],[],[],[],[]);
            w = (length(filter_weights{j})-1)/2;
            xobsj_sym = [flipud(xobs(2:w+1,j));xobs(:,j);flipud(xobs(end-w:end-1,j))];
            xobs(:,j) = conv(xobsj_sym,filter_weights{j},'valid');
        end

    elseif smoothing_window<0
        filter_weights = ones(-2*smoothing_window+1,1)/(-2*smoothing_window+1);
        xobs = [flipud(xobs(2:smoothing_window+1,:));xobs;flipud(xobs(end-smoothing_window:end-1,:))];
        xobs = conv2(filter_weights,1,xobs,'valid');
        filter_weights = repmat({filter_weights},1,n);
    else
        filter_weights = {};
    end

    %%% set integration line element
    dt = mean(diff(tobs));
    dv = mts*0+dt;%1./(2*mts+1);

    %%% get theta matrix
    Theta = build_theta(xobs,tags,custom_fcns,scale_x);

    %%% get column scaling 
    if scale_Theta == 0
        M_diag = ones(size(Theta,2),1);
    else
        M_diag = 1./prod(scale_x.^real(tags),2);
        cfun_offset = size(Theta,2)-length(M_diag);
        if cfun_offset
            M_diag = [M_diag;ones(cfun_offset,1)];
        end
    end

    grids = cell(n,1);
    Gs = cell(n,1);
    RTs = Gs;
    bs = Gs;
    vs = cell(n,2);
    bweaks = cell(n,1);

    for nn=1:n

        %%% get test function weights
        mt = mts(nn);
        pt = pts(nn);
        diffthresh = max(ceil((2*mt+1)*(1-overlap_frac)),1);
        [Cfs,~] = phi_int_weights(mt,max_d,-pt,phi_class);
        v = Cfs(end-1,:)*(mt*dt)^(-max_d+1)*dv(nn);
        vp = Cfs(end,:)*(mt*dt)^(-max_d)*dv(nn);
    
        %%% get test function centers
        [grid_i,bweak] = get_tf_centers(tobs,mt,K,diffthresh,relax_AG,xobs,nn,mts,pts,pt_ag_fac,mt_ag_fac,overlap_frac_ag);
        bweaks{nn} = bweak;

        %%% get linear system
        b = conv(xobs(:,nn),vp,'valid');
        b = b(grid_i);
        G = conv2(v,1,Theta,'valid');
        G = G(grid_i,:);

        %%% apply covariance
        if useGLS ~=0
            Vpmat = spdiags(repmat(vp/dv(nn),m-2*mt,1),0:2*mt,m-2*mt,m);
            Vpmat = Vpmat(grid_i,:);
            if useGLS<0
                Vmat = spdiags(repmat(v/dv(nn),m-2*mt,1),0:2*mt,m-2*mt,m);
                Vmat = Vmat(grid_i,:);
                alph_RT = -useGLS;
                Cov1 = Vpmat*Vpmat';
                Cov2 = Vmat*Vmat';
                Cov = (1-alph_RT)*Cov1 + alph_RT*Cov2;
            else
                alph_RT = useGLS;
                Cov = Vpmat*Vpmat';
    %             Cov = (1-alph_RT)*Cov + alph_RT*speye(size(Cov));
    %             Cov = Cov + alph_RT*speye(length(grid_i));
                Cov = (1-alph_RT)*Cov + alph_RT*diag(diag(Cov));
            end
            [RT,flag] = chol(Cov);
            if flag==0
                G = RT \ G;
                b = RT \ b;
            else
                disp('GLS failed: not pos. def')
            end
        else
            RT = [];
        end

        grids{nn} = grid_i;
        Gs{nn} = G;
        bs{nn} = b;
        RTs{nn} = RT;
        vs(nn,:) = {v,vp};
    end

end

function [grid_i,bweak] = get_tf_centers(tobs,mt,K,diffthresh,relax_AG,xobs,nn,mts,pts,pt_ag_fac,mt_ag_fac,overlap_frac_ag)

    m =length(tobs);
    dt = mean(diff(tobs));
    if relax_AG~=0
        if relax_AG>0
            diffthresh_ag = max(ceil((2*mt+1)*(1-overlap_frac_ag)),1);
            pt_ag = ceil(pts.*pt_ag_fac(1,:)+pt_ag_fac(2,:)); 
            mt_ag = ceil(mts.*mt_ag_fac(1,:)+mt_ag_fac(2,:)); 
            [Cfs_ag,~] = phi_int_weights(mt_ag(nn),1,-pt_ag(nn),1);        
            bweak = conv(xobs(:,nn),Cfs_ag(2,:)*(mt_ag(nn)*dt)^(-1),'valid');
            bweak = bweak(1:diffthresh_ag:end);
            Vmat = spdiags(repmat(Cfs_ag(1,:),m-2*mt_ag(nn),1),0:2*mt_ag(nn),m-2*mt_ag(nn),m);
            bweak = (Vmat(1:diffthresh_ag:end,:)')*((Vmat(1:diffthresh_ag:end,:)*Vmat(1:diffthresh_ag:end,:)')\bweak);
        else
            L = m-5*2;
            C = fdcoeffF(1, tobs(5+1), tobs(1:5*2+1));
            D = spdiags(repmat(C(:,end)',L,1), 0:5*2 ,L,m);
            bweak = D*xobs(:,nn);
            bweak = bweak(mts(nn)-5:end);
        end
        Y = abs(bweak);
        Y = cumsum(Y);
        Y = Y/Y(end);
        Y = abs(relax_AG)*Y+ (1-abs(relax_AG))*linspace(Y(1),Y(end),length(Y))';
        U = linspace(Y(1),Y(end),K);
%             U = sort(rand(1,K)*range(Y)+Y(1));
        grid_i = ones(1,K);
        for i = 2:K
            grid_temp = find(Y-U(i)>=0,1);
            df = grid_temp-grid_i(i-1);
            if df <= 0
                grid_i(i) = grid_i(i-1);
            elseif and(df < diffthresh,df>0)
                grid_i(i) = grid_i(i-1)+diffthresh;
            else
                grid_i(i) = grid_temp;
            end
        end
        grid_i = unique(min(max(grid_i-mt,1),m-2*mt));
    elseif relax_AG==0
        grid_i = 1:max(diffthresh,ceil((m-2*mt)/K)):m-2*mt;
        bweak = [];
    end

end
