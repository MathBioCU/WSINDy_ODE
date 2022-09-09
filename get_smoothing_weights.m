% W = get_smoothing_weights(xobs(:,2),tobs,m);
% 
% plot(W)

function W = get_smoothing_weights(xobs,tobs,m)

    dt = mean(diff(tobs));

    %%% get test function to estimate 2nd deriv
    [Cfs,~] = phi_int_weights(m,2,9,1);

    %%% estimate 2nd deriv in weak form
    fpp = rms(conv(xobs,Cfs(end,:)*(m*dt)^(-2)*dt,'valid'));

    %%% estimate variance
    sigma = estimate_sigma(xobs,mean(diff(tobs)));

    %%% 2nd moment weights
    b = fpp/2*dt^3*(-m:m)'.^2;

    cvx_solver SeDuMi
    cvx_precision([10^-16 10^-16 10^-16])
    cvx_begin quiet
        variable W(2*m+1,1)    
        minimize( sigma^2*dt^2*(W'*W) + b'*W )
        subject to
            {W >= 0, dt*sum(W) == 1};
    cvx_end

    W = W*dt;

end