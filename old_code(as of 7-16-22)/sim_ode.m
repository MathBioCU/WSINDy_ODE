%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% WSINDy: function for generating clean data
%%%%%%%%%%%% 
%%%%%%%%%%%% Copyright 2020, All Rights Reserved
%%%%%%%%%%%% Code by Daniel A. Messenger
%%%%%%%%%%%% For Paper, "Weak SINDy: Galerkin-based Data-Driven Model
%%%%%%%%%%%% Selection"
%%%%%%%%%%%% by D. A. Messenger and D. M. Bortz

function [weights,t,x,rhs] = sim_ode(x0,tspan,tol_ode,ode_name,params)

if strcmp(ode_name,'Logistic_Growth')
    pow = params{1};
    axi([2 pow+1],:) = [1;-1];
    rhs = @(x) x-x.^pow; 
    weights = {[[1 1];[pow -1]]};

elseif strcmp(ode_name,'Linear')
    A = params{1};
    axi = A';
    axi = [0*axi(1,:);axi];
    rhs = @(x) A*x; 
    weights = {};
    for dim = 1:size(A,1)
        weights{dim} = [eye(size(A,1)) A(dim,:)'];
    end
    
elseif strcmp(ode_name,'Duffing')
    mu = params{1}; alpha= params{2}; beta=params{3};
    axi(2:7,:) =[[0 1 0 0 0 0];[-alpha -mu 0 0 0 -beta]]';
    rhs = @(x) duff(x,mu,alpha,beta);
    weights = {[0 1 1], [[1 0 -alpha];[0 1 -mu];[3 0 -beta]]};

elseif strcmp(ode_name,'Lotka_Volterra')
    alpha = params{1}; beta = params{2}; delta = params{3}; gamma = params{4};
    axi(2:5,:) =[[alpha 0 0 -beta];[0 -gamma 0 delta]]';
    rhs = @(x) LoVo(x,alpha,beta,delta, gamma);
    weights = {[[1 0 alpha];[1 1 -beta]],[[0 1 -gamma];[1 1 delta]]};

elseif strcmp(ode_name,'Van_der_Pol')
    mu = params{1};
    axi(2:8,:) =[[0 1 0 0 0 0 0];[-1 mu 0 0 0 0 -mu]]';
    rhs = @(x) vanderpol(x,mu);
    weights = {[0 1 1], [[1 0 -1];[0 1 mu];[2 1 -mu]]};
    
elseif strcmp(ode_name,'Lorenz')
    sigma = params{1}; beta = params{2}; rho = params{3};
    axi(2:3,1) = [-sigma;sigma]; 
    axi(2:3,2) = [rho;-1]; 
    axi(4,3) = -beta;
    axi(7,2) = -1;
    axi(6,3) = 1;
    rhs = @(x) lorenz(x,sigma,beta,rho);    
    weights = {[[0 1 0 sigma];[1 0 0 -sigma]],...
                      [[1 0 0 rho];[1 0 1 -1];[0 1 0 -1]],...
                      [[1 1 0 1];[0 0 1 -beta]]};                  
elseif strcmp(ode_name,'Custom')
    rhs = params{1};
    weights = params{2};
end

options = odeset('RelTol',tol_ode,'AbsTol',tol_ode*ones(1,length(x0)));
[t,x]=ode45(@(t,x)rhs(x),tspan,x0,options);  

end


function dx = duff(x,mu,alpha,beta)
    dx(1) = x(2);
    dx(2) = -mu*x(2)-alpha*x(1)-beta*x(1).^3;
    dx = dx';
end

function dx = LoVo(x,alpha,beta,delta, gamma)
    dx(1) = alpha*x(1)-beta*x(1).*x(2);
    dx(2) = delta*x(1).*x(2)- gamma*x(2);
    dx = dx';
end

function dx = lorenz(x,sigma,beta,rho)
    dx(1) = sigma*(x(2)-x(1));
    dx(2) = x(1)*(rho-x(3))-x(2);
    dx(3) = x(1)*x(2)-beta*x(3);
    dx=dx';
end

function dx = vanderpol(x,mu)
    dx(1) = x(2);
    dx(2) = mu*x(2)-mu*x(1).^2.*x(2)-x(1);
    dx=dx';
end