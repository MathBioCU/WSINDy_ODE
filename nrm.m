function nrms = nrm(Ws,axi,p,inds)
    if isempty(inds)
	inds = 1:length(axi(:));
    end
    [J,N,M] = size(Ws);
    nrms = zeros(M,1);
    for m=1:M
	W = Ws(:,:,m);
	if p~=inf
            nrms(m) = norm(W(inds)-axi(inds),p)/norm(axi(inds),p);
	else
          nrms(m) = norm(abs(W(inds)-axi(inds))./abs(axi(inds)),p);
	end
    end
end

