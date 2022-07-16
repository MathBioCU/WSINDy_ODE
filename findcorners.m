function [mts,pts,sig_ests,corners] = findcorners(xobs,t,tau,tauhat,opt)

T = length(t);
n = size(xobs,2);
corners = zeros(n,2);
sig_ests = zeros(n,1);
mts = zeros(n,1);
pts = zeros(n,1);

for nn= 1:n
    [corner,sig_est] = findcornerpts(xobs(:,nn),t);
    k = corner(2);
    if opt ==1
        l = @(m,k,N) log((2*m-1)./m.^2).*(4*pi^2*k^2*m.^2-3*N^2*tauhat^2)-2*N^2*tauhat^2*log(tau);
        mnew = fzero(@(m)l(m,k,T), [1 2/sqrt(tau)]);
        if mnew>T/2-1
            mnew = T/2/k;
        end
        mts(nn) = min(floor((T-1)/2),ceil(mnew)); 
        pts(nn) = max(2,floor(log(tau)/log(1-(1-1/mts(nn))^2)));
    elseif opt == 2
        mnew = 1+T*tauhat/2/pi/k*sqrt(-2*log(tau));
        mts(nn) = min(floor((T-1)/2),ceil(mnew)); 
        pts(nn) = 2*pi*k/tauhat/T;
    end
    corners(nn,:)=corner;
    sig_ests(nn) = sig_est;
end
end

function [corner,sig_est] = findcornerpts(xobs,t)
    t = t(:);
    xobs = xobs(:);
    T = length(t);
    wn = ((0:T-1)-floor(T/2))'*(2*pi)/range(t);
    xx = wn(1:ceil(end/2));
    NN = length(xx); 
    Ufft = abs(fftshift(fft(xobs)));       
    Ufft = cumsum(Ufft);            
    Ufft = Ufft(1:ceil(T/2),:);
    errs = zeros(NN-2,1);
    for k=2:NN-1
       subinds1 = 1:k;
       subinds2 = k:NN;
       Ufft_av1 = Ufft(subinds1);
       Ufft_av2 = Ufft(subinds2);
       m1 = range(Ufft_av1)/range(xx(subinds1));
       m2 = range(Ufft_av2)/range(xx(subinds2)); 
       L1 = min(Ufft_av1)+m1*(xx(subinds1)-xx(1));
       L2 = max(Ufft_av2)+m2*(xx(subinds2)-xx(end));
       errs(k-1) = sqrt(sum(((L1-Ufft_av1)./Ufft_av1).^2) + sum(((L2-Ufft_av2)./Ufft_av2).^2));
    end
    [~,tstarind] = min(errs);
    tstar = -xx(tstarind);
    corner = [tstar max(NN-tstarind-3,1)];
    Ufft = fft(xobs)/sqrt(2*NN);
    sig_est = std(Ufft(max(floor(corner(2)*3),2):2*NN-floor(corner(2)*3)));
    sig_est = rms(Ufft(max(floor(corner(2)*3),2):2*NN-floor(corner(2)*3)));
end
