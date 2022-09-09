function [w_sparse,w_sparse_sindy,true_nz_weights,loss_wsindy,loss_sindy,ET_wsindy,ET_sindy,grids,pts,mts,Gs,bs,RTs,Theta_0,tags,bweaks,dxobs_0,vs] = ...
    wsindy_ode_fun(xobs,tobs,weights,...
    polys,trigs,...
    phi_class,max_d,tau,tauhat,K_frac,overlap_frac,relax_AG,...
    scale_Theta,useGLS,lambda,gamma,alpha_loss,...
    overlap_frac_ag,pt_ag_fac,mt_ag_fac,run_sindy,useFD,w)

    %%% get true weight vector
    n = size(xobs,2); 
    m = length(tobs);

    tags = get_tags(polys,trigs,n);
    J = size(tags,1);

    if ~isempty(weights)
        true_nz_weights = get_true_weights(weights,tags,n);
    else
        true_nz_weights = zeros(J,n);
    end

    %%% get scales
    if scale_Theta == 0
        M_diag = [];
        scale_x = ones(1,n);
    elseif scale_Theta < 0
        M_diag = ones(J,1);
        scale_x = ones(1,n);
    else
        scale_x = (vecnorm(xobs.^max(real(tags)),2)).^(1./max(real(tags)));    
        M_diag = 1./prod(scale_x.^real(tags),2);        
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
            [mts,pts,sig_ests,corners] = findcorners(xobs,tobs,tau,-tauhat,phi_class);
        else
            disp('Error: tauhat cannot be 0')
            return
        end
    end
    
    %%% set integration line element
    dt = mean(diff(tobs));
    dv = mts*0+dt;%1./(2*mts+1);

    %%% smooth data
    if w>0
        smoothing_grid = linspace(-1,1,2*w+3);
        filter_weights = 1-smoothing_grid.^2;
        filter_weights = filter_weights/sum(filter_weights);
        filter_weights = filter_weights(2:end-1);
        xobs = [flipud(xobs(2:w+1,:));xobs;flipud(xobs(end-w:end-1,:))];
        xobs = conv2(filter_weights,1,xobs,'valid');
    elseif w<0
        xobs = movmean(xobs,-w,1,'Endpoints','shrink');
    end

    %%% get theta matrix
    Theta_0 = build_theta(xobs,tags,scale_x);

    w_sparse = zeros(size(Theta_0,2),n);
    grids = cell(n,1);
    Gs = cell(n,1);
    RTs = Gs;
    bs = Gs;
    vs = cell(n,2);
    bweaks = cell(n,1);
    loss_wsindy = cell(n,1);

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
        b = conv(xobs(:,nn)/scale_x(nn),vp,'valid');
        b = b(grid_i);
        G = conv2(v,1,Theta_0,'valid');
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

        %%% regress
        if length(lambda)==1
            [w_sparse(:,nn),its(nn)] = sparsifyDynamics(G,b,lambda,1,gamma,M_diag*scale_x(nn));
            loss_wsindy = [];
        else
            if ~isempty(M_diag)
                M_scale_b = [1/scale_x(nn);M_diag];
            else
                M_scale_b = [];
            end
            [w_sparse(:,nn),loss_wsindy{nn},its(nn)] = wsindy_pde_RGLS_seq(lambda,gamma,[b G],1,M_scale_b,alpha_loss);
        end

        %%% save outputs
        grids{nn} = grid_i;
        Gs{nn} = G;
        bs{nn} = b;
        RTs{nn} = RT;
        vs(nn,:) = {v,vp};
    end
    ET_wsindy = toc;

    %%% -----------------SINDy--------------------------- 

    tic;
    if run_sindy
    loss_sindy = cell(n,1);
    
    if useFD<=0
        L = length(tobs)-2;
        C = fdcoeffF(1, tobs(2), tobs(1:3));
        D = spdiags(repmat(C(:,end)',L,1), 0:2 ,L,length(tobs));
        naive_diff = D*xobs;
        dxobs_0 = 0*naive_diff;
        if useFD ==0
            reg_param = 1/sqrt(length(tobs));
        else
            reg_param = (-useFD)/sqrt(length(tobs));
        end
        for j=1:n
            dxobs_0(:,j) = TVRegDiff( xobs(2:end-2,j), 20, reg_param, naive_diff(:,j), [], [], dt, 0, 0 );
        end
        Theta_sindy = Theta_0(max(useFD,1)+1:end-max(useFD,1),:);
    else
        L = length(tobs)-useFD*2;
        C = fdcoeffF(1, tobs(useFD+1), tobs(1:useFD*2+1));
        D = spdiags(repmat(C(:,end)',L,1), 0:useFD*2 ,L,length(tobs));
        Cov_s = D*D';
        useGLS = 0; %%% uncomment to apply covariance to standard sindy -  results improve!
        if useGLS~=0
            alph_RT = abs(useGLS);
            Cov_s = (1-alph_RT)*Cov_s + alph_RT*diag(diag(Cov_s));
            [RT_s,flag] = chol(Cov_s);
            dxobs_0 = RT_s \ (D*xobs);
            Theta_sindy = RT_s \ Theta_0(max(useFD,1)+1:end-max(useFD,1),:);   
        else
            dxobs_0 = D*xobs;
            Theta_sindy = Theta_0(max(useFD,1)+1:end-max(useFD,1),:);   
        end
    end

    w_sparse_sindy = 0*w_sparse;

    for nn=1:n
        if length(lambda)==1
            [w_sparse_sindy(:,nn),its(nn)] = sparsifyDynamics(Theta_sindy,dxobs_0(:,nn),lambda,1,gamma,M_diag(2:end));
            loss_sindy = [];
        else            
            if ~isempty(M_diag)
                M_scale_b = [1;M_diag];
            else
                M_scale_b = [];
            end
            [w_sparse_sindy(:,nn),loss_sindy{nn},its(nn)] = wsindy_pde_RGLS_seq(lambda,gamma,[dxobs_0(:,nn) Theta_sindy],1,M_scale_b,alpha_loss);
        end
    end
    ET_sindy = toc;
    else
        w_sparse_sindy = 0*true_nz_weights;
        ET_sindy = 0;
        loss_sindy = loss_wsindy;
        dxobs_0 = [];
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