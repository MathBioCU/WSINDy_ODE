function Tps = tpscore(Ws,axi)
    [J,N,M] = size(Ws);
    tnz = find(axi);
    Tps = zeros(M,1);
    for m=1:M
	nz = find(Ws(:,:,m));
	FN = length(setdiff(tnz,nz));
	FP = length(setdiff(nz,tnz));
	TP = length(intersect(tnz,nz));
	Tps(m) = TP/(TP+FN+FP);
    end
end
