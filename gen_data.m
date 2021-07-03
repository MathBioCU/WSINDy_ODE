function [weights,x,t,x0,ode_name,ode_params] = gen_data(ode_num,ode_params,tspan,x0,tol_ode)

% Simulate ODE
    ode_names = {'Linear','Logistic_Growth','Van_der_Pol','Duffing','Lotka_Volterra','Lorenz','Rossler','rational','Oregonator','Hindmarsh-Rose','Pendulum','custom'};
    ode_name = ode_names{ode_num};
    if strcmp(ode_name, 'Linear')%1
        if isempty(ode_params)
            ll = 5;
            x = zeros(ll,1);y= zeros(ll,1); x([2 end])=[1 -1]; y([2 end]) = [-1 1];
            A = toeplitz(y,x');
            A(logical(eye(size(A,1))))= -0.2;
            ode_params = {A};
        end
        if isempty(x0)
            x0 = zeros(ll,1); 
            x0(1) = 10; 
        end
        if isempty(tspan)
            tspan = 0:0.025:40;
        end

    elseif strcmp(ode_name, 'Logistic_Growth')%2
        if isempty(ode_params)
            ode_params = {2};
        end
        if isempty(x0)
            x0=0.01;
        end
        if isempty(tspan)
            tspan =  0:0.005:10;
        end

    elseif strcmp(ode_name, 'Van_der_Pol')%3
        if isempty(ode_params)
            ode_params = {4};
        end
        if isempty(x0)
            x0=[0;1];
        end
        if isempty(tspan)
            tspan = 0:0.01:30;    
        end

    elseif strcmp(ode_name, 'Duffing')%4
        if isempty(ode_params)
            mu = 0.2; ode_params = {mu, mu, 1}; 
        end
        if isempty(x0)
            x0=[0;2];
        end
        if isempty(tspan)
            tspan = 0:0.01:30;    
        end
        
    elseif strcmp(ode_name, 'Lotka_Volterra')%5
        if isempty(ode_params)
            alpha = 3; beta = 1; ode_params = {alpha,beta,beta,2*alpha};     
        end
        if isempty(x0)
            x0=[1;1];
        end
        if isempty(tspan)
            tspan = 0:0.01:10;
        end
        
    elseif strcmp(ode_name, 'Lorenz')%6
        %x0 = [5 5 25]';%x0 = [rand(2,1)*30-15;rand*30+10]
        if isempty(ode_params)
            ode_params = {10, 8/3, 28};     
        end
        if isempty(x0)
            x0=[-8 7 27]';
        end
        if isempty(tspan)
            tspan = 0:0.001:10;
        end

    elseif strcmp(ode_name, 'Rossler')%7
    %        ode_params = {0.2, 0.2, 5.7}; x0 = [3 5 0]'; tspan = 0:.005:25;
        if isempty(ode_params)
            ode_params = {0.2, 0.2, 5.7}; 
        end
        if isempty(x0)
            x0 = [3 5 0]';
        end
        if isempty(tspan)
            tspan = 0:.01:25;
        end
        
    elseif strcmp(ode_name, 'rational')%8
        if isempty(ode_params)
            ode_params = {8,1}; 
        end
        if isempty(x0)
            x0 = -5;
        end
        if isempty(tspan)
            tspan = 0:.005:20;
        end
        

    elseif strcmp(ode_name, 'Oregonator')%9
        if isempty(ode_params)
            ode_params = {0.5, 5, 5, 0.5, 0.5, 0.95};
        end
        if isempty(x0)
            x0 = [1 0.2 50]';    
        end
        if isempty(tspan)
            tspan = 0:.002:50;
        end
                
    elseif strcmp(ode_name, 'Hindmarsh-Rose')%10
        if isempty(ode_params)
            a = 1; b = 3; c = 1; 
            d = 5; r = 10^-3; s = 4;
            xR = -3.19/4; I = 0;
            ts = 10;
            ode_params = {a,b,c,d,r,s,xR,I,ts};
        end
        if isempty(x0)
            x0 = [-1.3095   -7.5901   -0.2020]';
        end
        if isempty(tspan)
            tspan = 0:.001:100;
        end

    elseif strcmp(ode_name, 'Pendulum')%11
        if isempty(ode_params)
            ode_params = {1};
        end
        if isempty(x0)
            x0 = [pi-1*pi/16 0]';
        end
        if isempty(tspan)
            tspan = 0:0.1:50;
        end

    elseif strcmp(ode_name, 'custom')%12
        % nonlinear Schrodinger
    %     ode_func = @(x) [0.1*x(2).*(x(1).^2+x(2).^2);-x(1).*(x(1).^2+x(2).^2)];
    %     weights = {[[2 1 0.1];[0 3 0.1]], [[3 0 -1];[1 2 -1]]};
        % cubic oscillator
        if isempty(ode_params)
            ode_func = @(x) [-0.1*x(1).^3+2*x(2).^3;-2*x(1).^3-0.1*x(2).^3];
            weights = {[[3 0 -0.1];[0 3 2]], [[3 0 -2];[0 3 -0.1]]};
            ode_params = {ode_func,weights};
        end
        if isempty(x0)
            x0 = [1;0];
        end
        if isempty(tspan)
            tspan = 0:.01:20;
        end        

    end
    [weights,t,x,~] = sim_ode(x0,tspan,tol_ode,ode_name,ode_params,'o');

end

