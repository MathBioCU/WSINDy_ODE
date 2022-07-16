function true_nz_weights = get_true_weights(weights,tags,n)
    true_nz_weights = zeros(size(tags,1),n);
    for i = 1:length(weights)
        weights_i = weights{i};
        [l1,l2] = size(weights_i);
        for j = 1:l1
            true_nz_weights(all(weights_i(j,1:l2-1) == tags,2),i) = weights_i(j,l2);
	end
    end
end
