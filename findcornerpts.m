function [corner,sig_est] = findcornerpts(xobs,t)
    t = t(:);
    T = length(t);
    wn = ((0:T-1)-floor(T/2))'*(2*pi)/range(t);
    xx = wn(1:ceil(end/2));
    NN = length(xx); 
    Ufft = mean(abs(fftshift(fft(xobs))),2)/sqrt(2*NN);
    Ufft = Ufft(1:ceil(T/2),:);
    [~,Umax] = max(Ufft);
    tstarind1 = getcorner(cumsum(abs(Ufft(1:Umax))),xx(1:Umax));
    tstarind2 = getcorner(log(abs(Ufft(1:Umax))),xx(1:Umax));
    tstarind = floor((tstarind1+tstarind2)/2);
    %%% uncomment to use matlab's changepoint algorithm on either the
    %%% cumulative sum or the log of the power spectrum:
    %     tstarind = findchangepts(cumsum(Ufft));
    %     tstarind = findchangepts(log(Ufft));
    tstar = -xx(tstarind);
    corner = [tstar max(Umax-tstarind+1,1)];
    sig_est = rms(Ufft(1:min(corner(:,2)),:));
end

function tstarind = getcorner(Ufft,xx)
    NN = length(Ufft);
    Ufft = Ufft/max(abs(Ufft))*NN;
    errs = zeros(NN,1);
    for k=1:NN
        [L1,L2,m1,m2,b1,b2,Ufft_av1,Ufft_av2]=build_lines(Ufft,xx,k);
%         plot(1:k,L1,k:length(xx),L2,1:length(xx),Ufft,'b-o'); drawnow; % uncomment to visualize alg
        errs(k) = sqrt(sum(((L1-Ufft_av1)./Ufft_av1).^2) + sum(((L2-Ufft_av2)./Ufft_av2).^2)); % pointwise relative 2-norm 
%         errs(k) = m1^2/(1+m1^2)*norm(1:k)^2 + 1/(1+m1^2)*norm(Ufft_av1-b1)^2 + m2^2/(1+m2^2)*norm(k:NN)^2 + 1/(1+m2^2)*norm(Ufft_av2-b2)^2; % total least squares v1
%         errs(k) = norm(L1-Ufft_av1,2)/norm(Ufft_av1,2)+norm(L2-Ufft_av2,2)/norm(Ufft_av2,2); % relative 2-norm 
%         errs(k) = norm(L1-Ufft_av1,inf)/norm(Ufft_av1,inf)+norm(L2-Ufft_av2,inf)/norm(Ufft_av2,inf); % relative inf-norm 
%         A = [ones(NN+1,1) [-m1*ones(length(Ufft_av1),1);-m2*ones(length(Ufft_av2),1)]]; % total least squares v2 
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