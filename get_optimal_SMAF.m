% get optimal simple moving average filter, based on mean-squared error minimization
function [W,sigma_est,its] = get_optimal_SMAF(x,fx_obs,max_points,init_m_fac,max_filter_fac,expand_fac,maxits,deriv_tol,verbose)

    if isempty(max_points)
        max_points = 10^5;
    end
    if isempty(init_m_fac)
        init_m_fac = 200;
    end
    if isempty(max_filter_fac)
        max_filter_fac = 8;
    end
    if isempty(expand_fac)
        expand_fac = 2;
    end
    if isempty(maxits)
        maxits = 100;
    end
    if isempty(deriv_tol)
        deriv_tol = 10^-6;
    end
    if isempty(verbose)
        verbose = 0;
    end

    sigma_est = estimate_sigma(fx_obs);
    subsamp = max(floor(length(x)/max_points),1);  
    
    fx_subsamp = fx_obs(1:subsamp:end,:);
    M = size(fx_subsamp,1);
    dx_subsamp = mean(diff(x(1:subsamp:end)));
    
    m = ceil(M/init_m_fac);
    max_filter_width=floor(M/max_filter_fac);
    
    its = 1; check=1; 
    m = min(m,max_filter_width);

    while and(check>0,its<maxits)
        if verbose
            tic
        end
        [~,A] = build_poly_kernel(2,@(x) x*0+1,min(max(floor(m*expand_fac),3),floor((M-1)/2)),dx_subsamp,0);
        if size(fx_subsamp,2)==1
            d = 2*mean(abs(conv(fx_subsamp,A(3,:),'valid')));
        else
            d = 2*mean(reshape(abs(conv2(A(3,:),1,fx_subsamp,'valid')),[],1));
        end
        C = sigma_est^2/((d+deriv_tol)^2*dx_subsamp^4/144);
        mnew = min( floor((fzero(@(n) n.^5-n.^3-C,1)-1)/2), max_filter_width);
    
        check = abs(m-mnew);
        m = mnew;
        its = its+1;
        if verbose
            disp([toc m d])
        end
    end
    m = m*subsamp;
    W = 1./(2*m+1)*ones(2*m+1,1);
    
end

function [f,A]=build_poly_kernel(deg,k,n,dx,max_dx)
    x = (-n:n)'*dx;
    X = x.^(0:deg);
    K = k(x/(n*dx));
    K = K/norm(K,1);
    A = pinv(sqrt(K).*X).*sqrt(K)';
    M = [diag(factorial(0:max_dx)) zeros(max_dx+1,deg-max_dx)];
    f = M*A;
end