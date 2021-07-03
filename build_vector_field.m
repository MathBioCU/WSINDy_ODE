function rhs_dd = build_vector_field(coeffs,tags)
    
dim = size(coeffs,2);
rhs = cell(dim,1);

for j=1:dim
    inds = find(coeffs(:,j));
    g_temp = @(x) 0;
    for i=1:length(inds)
        if isreal(tags(inds(i),:))
            g_temp = @(x) g_temp(x) + coeffs(inds(i),j)*prod(x.^(tags(inds(i),:)'),1);
        else
            imag_tags = imag(tags(inds(i),:));
            ind_temp = find(imag_tags);
            if imag_tags(ind_temp)<0
                 g_temp = @(x) g_temp(x) + coeffs(inds(i),j)*sin(x(ind_temp)*abs(imag_tags(ind_temp)));
            else
                 g_temp = @(x) g_temp(x) + coeffs(inds(i),j)*cos(x(ind_temp)*abs(imag_tags(ind_temp)));
            end
       end
    end
    rhs{j} = g_temp;
end

rhs_dd = @(x) rhs{1}(x);
for k=2:dim 
    rhs_dd = @(x) [rhs_dd(x);rhs{k}(x)];
end

end    
