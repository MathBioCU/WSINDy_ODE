function Tps = tpscore(Ws,axi)
    tnz = find(axi);
    if isequal(class(Ws),'cell')
        M=length(Ws);
        Tps = zeros(M,1);
        for m=1:M
            nz = find(Ws{m});
            FN = length(setdiff(tnz,nz));
            FP = length(setdiff(nz,tnz));
            TP = length(intersect(tnz,nz));
            Tps(m) = TP/(TP+FN+FP);
        end
    elseif isequal(class(Ws),'double')
        dims_axi = size(axi);
        M=numel(Ws)/numel(axi);
        Tps = zeros(M,1);
        for m=1:M
            if dims_axi(2)>1
                nz = find(Ws(:,:,m));
            else
                nz = find(Ws(:,m));
            end
            FN = length(setdiff(tnz,nz));
            FP = length(setdiff(nz,tnz));
            TP = length(intersect(tnz,nz));
            Tps(m) = TP/(TP+FN+FP);
        end
    else
        Tps=NaN;
    end
end
