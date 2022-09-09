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