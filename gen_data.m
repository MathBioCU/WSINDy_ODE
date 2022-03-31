if strcmp(ode_name, 'Linear')
    if use_preset_params
        % set params common to WSINDy and SINDy 
        polys = 0:5;                    % monomial powers to include in library
        trigs = [];                     % sine / cosine frequencies to include in library
        gamma = 10^-Inf;                % Tikhonoff regularization parameter
        lambda_mult = 4;                % sets sparsity knob lambda = min(true_weights)/lambda_mult; lambda_mult = Inf => lambda = 0
        scale_Theta = 0;                % toggle normalize columns of Theta_0 by norm = scale_Theta_0

        % set WSINDy params            
        tau = 1;                        % toggle adaptive grid - convex comb parameter between uniform (=0) and adaptive (=1)
        tau_p = 16;                     % test function has value 10^-tau_p at penultimate support point. Or, if tau_p<0, directly sets poly degree p = -tau_p
        K = 126;                        % num basis fcns
        p = 2; s = 16;                  % AG weak-deriv params -- poly degree, support size
        r_whm = 30;                     % width at half-max in num of timepoints (s.t. phi(r_whm*dt) = 1/2) 
        useGLS = 10^-12;                % toggle use generalized least squares
    end
elseif strcmp(ode_name, 'Logistic_Growth')
    if use_preset_params
        % set common params
        polys = 0:5;                    
        trigs = [];                    
        gamma = 10^-2.0;                
        lambda_mult = 4;                
        scale_Theta = 2;               

        % set WSINDy params             
        tau = 1;                        
        tau_p = 16;                     % test function has value 10^-tau_p at penultimate support point. Or, if tau_p<0, directly sets poly degree p = -tau_p
        K = 100;                        
        p = 2; s = 16;                  
        r_whm = 30;                     
        useGLS = 10^-12;                % toggle use generalized least squares
    end
elseif strcmp(ode_name, 'Van_der_Pol')
    if use_preset_params
        % set common params
        polys = 0:5;                   
        trigs = [];                     
        gamma = 10^-Inf;
        lambda_mult = 4;
        scale_Theta = 0;

        % set WSINDy params
        tau = 1;
        tau_p = 16;                     % test function has value 10^-tau_p at penultimate support point. Or, if tau_p<0, directly sets poly degree p = -tau_p
        K = 126;
        p = 2; s = 16;
        r_whm = 30; 
        useGLS = 10^-12;                % toggle use generalized least squares
    end
elseif strcmp(ode_name, 'Duffing')
    if use_preset_params
        % set common params
        polys = 0:5;
        trigs = [];
        gamma = 10^-Inf;
        lambda_mult = 2;
        scale_Theta = 0;

        % set WSINDy params
        tau = 1;
        tau_p = 16;                     % test function has value 10^-tau_p at penultimate support point. Or, if tau_p<0, directly sets poly degree p = -tau_p
        K = 126;
        p = 2; s = 16;
        r_whm = 30;
        useGLS = 10^-12;                % toggle use generalized least squares
    end
elseif strcmp(ode_name, 'Lotka_Volterra')
    if use_preset_params
        % set common params
        polys = 0:5;
        trigs = [];
        gamma = 0;
        lambda_mult = 4;
        scale_Theta = 0;

        % set WSINDy params
        tau = 1;
        tau_p = 16;                     % test function has value 10^-tau_p at penultimate support point. Or, if tau_p<0, directly sets poly degree p = -tau_p
        K = 126;
        p = 2; s = 16;
        r_whm = 30;
        useGLS = 10^-12;                % toggle use generalized least squares
    end
elseif strcmp(ode_name, 'Lorenz')
    if use_preset_params
        % set common params
        polys = 0:5;
        trigs = [];
        gamma = 10^-Inf;
        lambda_mult = 8;
        scale_Theta = 0;

        % set WSINDy params
        tau = 1;
        tau_p = 16;                     % test function has value 10^-tau_p at penultimate support point. Or, if tau_p<0, directly sets poly degree p = -tau_p
        K = 224;
        p = 2; s = 16;
        r_whm = 30;
        useGLS = 10^-12;                % toggle use generalized least squares
    end
end


n = length(x0); 
[weights,t,x,rhs] = sim_ode(x0,tspan,tol_ode,ode_name,ode_params);

%%% add noise

signal_power = rms(x(:));
sigma = noise_ratio*signal_power;
noise = normrnd(0,sigma,size(x));
xobs = x + noise;
tobs = t;

