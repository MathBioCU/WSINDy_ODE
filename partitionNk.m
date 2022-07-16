function total_parts = partitionNk(N,k)
    if N>1
        A = partitions(N,0:N,N,k);
	num_unique_parts = size(A,1);
        unique_parts = zeros(num_unique_parts,k);
	for j=1:num_unique_parts
            temp_part = [];
            for m= 1:N+1
		if A(j,m)>0
                    for l = 1:A(j,m)
                        temp_part = [temp_part m-1];
                    end
                end
            end
            unique_parts(j,:) = temp_part;
        end
        total_parts = [];
        for j=1:num_unique_parts
            total_parts = [total_parts;unique(perms(unique_parts(j,:)),'rows')];
        end
    elseif N==1
        total_parts = eye(k);
    elseif N==0
        total_parts = zeros(1,k);
    end
end
