function Cfs = phi_int_weights(mt,pt)
    %%% Cfs are the coefficients of the test function 
    %%% phi(t) = (1-(t/(mt*dt))^2)^(pt)
    %%% and its derivatives, which can be computed manually but when 
    %%% the differentiation order is large its easier to use this code
    
    d = 1; % change for higher-order systems
    t = (0:mt)/mt;
    t_L = zeros(d+1,mt+1);                             % store (1+t)^q, (1-t)^q
    t_R = zeros(d+1,mt+1);                                  
    for j=1:mt
        t_L(:,j)  = (1+t(j)).^(fliplr(pt-d:pt))';         
        t_R(:,j)  = (1-t(j)).^(fliplr(pt-d:pt))';
    end
    ps = ones(d+1,1);                                  % derivative coefficients
    for q=1:d
        ps(q+1) = (pt-q+1)*ps(q);
    end
    t_L = ps.*t_L;
    t_R = ((-1).^(0:d)'.*ps).*t_R;

    Cfs = zeros(d+1,2*mt+1);                            % Values of derivatives at grid points
    Cfs(1,:) = [fliplr(t_L(1,:).*t_R(1,:)) t_L(1,2:end).*t_R(1,2:end)];
    P = fliplr(pascal(d+1));    
    for k=1:d
        binoms = diag(P,d-k);
        Cfs_temp = zeros(1,mt+1);
        for j=1:k+1
            Cfs_temp = Cfs_temp + binoms(j)*t_L(k+2-j,:).*t_R(j,:);
        end
        Cfs(k+1,:) = [(-1)^k*fliplr(Cfs_temp) Cfs_temp(2:end)];
    end
end    
