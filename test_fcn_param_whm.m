function [p,m] = test_fcn_param_whm(r,t,tol)
    dt = t(2)-t(1);
    r_whm = r*dt;
    A = log2(10)*(-log10(tol));

    gg = @(s) -s.^2.*((1-(r_whm./s).^2).^A-1);
    hh = @(s) (s-dt).^2;
    ff = @(s) hh(s)-gg(s);

    s = fzero(ff,[r_whm, r_whm*sqrt(A)+dt]);
    m = ceil(s/dt);
    p = min(ceil(max(-1/log2(1-(r_whm/s)^2),1)),200);
end