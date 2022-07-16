common_params = {polys,trigs,lambda_mult,scale_Theta,gamma};
wsindy_params = {s, K, p, tau};
[Theta_0, tags, true_nz_weights, M_diag, lambda] = build_theta(xobs,weights,common_params,n);

tic;
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
    [V,Vp,ab_grid,ps] = VVp_build_adaptive_whm(tobs,grid_i, r_whm, tau_p, {0,inf,0});  %{pow,nrm,ord}. ord=0, ||phi||, ord=1, ||phi'|| 
    ps_all = [ps_all;ps];
    mats{i} = {V,Vp};
    ts_grids{i} = ab_grid;
    Ys{i} = Y;    
    if useGLS > 0
        Cov = Vp*Vp'+useGLS*eye(size(V,1));
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
