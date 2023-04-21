function [U,noise,noise_ratio_obs,sigma] = gen_noise(U_exact,sigma_NR,noise_dist,noise_alg)

    if noise_alg == 0 % additive
        stdv = rms(U_exact(:)).^2;
    elseif noise_alg == 1 % multiplicative
        stdv = 1;
    end
    dims = size(U_exact);
    if noise_dist == 0 % white noise
        if sigma_NR>0
            sigma = sigma_NR*sqrt(stdv);
        else
            sigma=-sigma_NR;
        end
        noise = normrnd(0,sigma,dims);
    elseif noise_dist == 1 % uniform noise
        if sigma_NR>0
            sigma = (3*sigma_NR^2*stdv)^(1/2);
        else
            sigma=-sigma_NR;
        end
        noise = sigma*(2*rand(dims)-1);
    end
    if noise_alg == 0 % additive
        U = U_exact + noise;
    elseif noise_alg == 1 % multiplicative
        U = U_exact.*(1 + noise);
    end
    noise_ratio_obs = norm(U(:)-U_exact(:))/norm(U_exact(:));
    
end