function [corner,sig_est] = findcornerpts(xobs,t)
    t = t(:);
    xobs = xobs(:);
    T = length(t);
    wn = ((0:T-1)-floor(T/2))'*(2*pi)/range(t);
    xx = wn(1:ceil(end/2));
    NN = length(xx); 
    Ufft = fftshift(fft(xobs))/sqrt(2*NN);
    tstarind1 = getcorner(cumsum(abs(Ufft(1:ceil(T/2),:))),xx);
    tstarind2 = getcorner(log(abs(Ufft(1:ceil(T/2),:))),xx);
    tstarind = floor((tstarind1+tstarind2)/2);
    tstar = -xx(tstarind);
    corner = [tstar max(NN-tstarind-3,1)];
    Ufft=fftshift(Ufft);
    sig_est = rms(Ufft(max(floor(corner(2)*3),2):2*NN-floor(corner(2)*3)));
end

function tstarind = getcorner(Ufft,xx)
    NN = length(Ufft);
    Ufft = Ufft/max(abs(Ufft))*NN;
    errs = zeros(NN,1);
    for k=1:NN
        [L1,L2,m1,m2,b1,b2,Ufft_av1,Ufft_av2]=build_lines(Ufft,xx,k);
        errs(k) = sqrt(sum(((L1-Ufft_av1)./Ufft_av1).^2) + sum(((L2-Ufft_av2)./Ufft_av2).^2)); % relative l2 
%             errs(k) = m1^2/(1+m1^2)*norm(1:k)^2 + 1/(1+m1^2)*norm(Ufft_av1-b1)^2 + m2^2/(1+m2^2)*norm(k:NN)^2 + 1/(1+m2^2)*norm(Ufft_av2-b2)^2; % relative l2 
%            errs(k) = norm(L1-Ufft_av1,2)/norm(Ufft_av1,2)+norm(L2-Ufft_av2,2)/norm(Ufft_av2,2); % relative l2 
%            errs(k) = norm(L1-Ufft_av1,inf)/norm(Ufft_av1,inf)+norm(L2-Ufft_av2,inf)/norm(Ufft_av2,inf); % relative l2 
%         A = [ones(NN+1,1) [-m1*ones(length(Ufft_av1),1);-m2*ones(length(Ufft_av2),1)]];
%         c = [-b1*ones(length(Ufft_av1),1);-b2*ones(length(Ufft_av2),1)];
%         z = [[Ufft_av1(:) (1:k)']*(m1-1)/(m1^2+1); [Ufft_av2(:) (k:NN)']*(m2-1)/(m2^2+1)];
%         errs(k) = norm(sum(A.*z,2)-c);        
        
    end
    [~,tstarind] = min(errs);
end

function [L1,L2,m1,m2,b1,b2,Ufft_av1,Ufft_av2]=build_lines(Ufft,xx,k)

   NN = length(Ufft);
   subinds1 = 1:k;
   subinds2 = k:NN;
   Ufft_av1 = Ufft(subinds1);
   Ufft_av2 = Ufft(subinds2);
   
   [m1,b1,L1]=lin_regress(Ufft_av1,xx(subinds1));
   [m2,b2,L2]=lin_regress(Ufft_av2,xx(subinds2));

end

function [m,b,L]=lin_regress(U,x)

   m = (U(end)-U(1))/(x(end)-x(1));
   b = U(1)-m*x(1);
   L = U(1)+m*(x-x(1));

%    A = [0*x(:)+1 x(:)];
%    c = A \ U(:);
%    b = c(1); m = c(2); L = A*c;

end