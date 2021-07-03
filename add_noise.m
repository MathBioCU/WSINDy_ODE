
function [sigma,noise,xobs,tobs,noise_ratio_obs] = add_noise(x,noise_ratio,t)
    signal_power = rms(x(:));
    if noise_ratio<0
        sigma = -noise_ratio;
    else
        sigma = noise_ratio*signal_power;
    end
    noise = normrnd(0,sigma,size(x));
    xobs = x + noise;
    tobs = t;
    noise_ratio_obs = norm(noise(:))/norm(x(:));
end