w = w_sparse;
w = w_sparse_sindy;
fcns = custom_fcns;

w_params = cell(1,size(w,2));
w_features = cell(size(w,2),1);
for i=1:size(w,2)
    w_params{i} = w(w(:,i)~=0,i);
    w_features{i} = fcns(w(:,i)~=0);
end
x0_test = x(1,:);
t_test = tobs;

rhs_p = @(x,params) rhs_fun(w_features,params,x);

tol_ode = 1e-12;
options_ode_sim = odeset('RelTol',tol_ode,'AbsTol',tol_ode*ones(1,length(x0_test)));
[t_learned,x_learned]=ode45(@(t,x)rhs_p(x,w_params),tobs,x0_test,options_ode_sim);

plot(t_learned,x_learned,tobs,xsub)
title(norm(x_learned(:)-xsub(:))/norm(xsub(:)))

function dx = rhs_fun(features,params,x)
    nstates = length(x);
    x = num2cell(x);
    dx = zeros(nstates,1);
    for i=1:nstates
        dx(i) = cellfun(@(z1) z1(x{:}),features{i})*params{i}(:);
    end
end
