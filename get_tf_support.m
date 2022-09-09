function m = get_tf_support(phi,N,tauhat,k_x)
    ks = ((0:N-1)-floor(N/2)).^2;
    errs = zeros(floor(N/2)-1,1);
    m = 1;
    errs(1) = abs((k_x/tauhat)^2-sum(ks));
    check = 0;
    while and(check == 0,m<=(N-3)/2)
        m = m + 1;
        x = -1:(1/m):1;
        phi_grid = [phi(x) zeros(1,N-2*m-1)];
        phi_fft = fftshift(abs(fft(phi_grid)));
        phi_fft = phi_fft/sum(phi_fft);
        errs(m) = abs((k_x/tauhat)^2-sum(phi_fft.*ks));
        check = errs(m)>errs(m-1);
    end
end