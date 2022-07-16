function Theta = build_theta(xobs,tags,scale_x)

[m,n] = size(xobs);
J = size(tags,1);
Theta = zeros(m,J);
if isempty(scale_x)
    scale_x = ones(1,n);
end

for j = 1:J
    if all(isreal(tags(j,:)))
        Theta(:,j) = prod((xobs./scale_x).^(tags(j,:)),2);
    else
        tag = tags(j,:);
        ind = find(tag);
        if imag(tag(ind))<0
            Theta(:,j) = sin(-imag(tag(ind))*xobs(:,ind));
        else
            Theta(:,j) = cos(imag(tag(ind))*xobs(:,ind));
        end
    end
end
end