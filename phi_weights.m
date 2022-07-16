function Cfs = phi_weights(phifun,m,maxd)
    x = linspace(-1,1,2*m+1);
    Cfs = zeros(maxd+1,2*m+1);
    syms y;
    f = @(y)phifun(y);
    for j=1:maxd+1
        Df = matlabFunction(diff(f(y),j-1));
        Cfs(j,:) = fillmissing(Df(x),'constant',Df(eps));
        inds = find(isinf(abs(Cfs(j,:))));
        for k=1:length(inds)
            Cfs(j,inds(k)) = Df(x(inds(k))+eps);
        end
    end
    Cfs = Cfs/norm(Cfs(1,:),1);
end
