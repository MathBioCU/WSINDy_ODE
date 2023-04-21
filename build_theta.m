function Theta = build_theta(xobs,tags,custom_fcns,scale_x)

    [M,nstates] = size(xobs);
    if isempty(scale_x)
        scale_x = ones(1,nstates);
    end
    Theta_0 = arrayfunvec(tags',@(tag)get_col(xobs,tag,scale_x),1);

    xobs_cell = mat2cell(xobs,M,ones(1,nstates));
    try
        Theta_1 = cell2mat(cellfun(@(y) y(xobs_cell{:}), custom_fcns(:), 'uni',0)');
    catch
        disp('custom fcs improperly specified, omitting from library')
        Theta_1 = [];
    end
    Theta = [Theta_0 Theta_1];
 
end

function col = get_col(xobs,tag,scale_x)
    tag=tag(:)';
    if all(isreal(tag))
        col = prod((xobs./scale_x).^tag,2);
    else
        ind = find(tag);
        if imag(tag(ind))<0
            col = sin(-imag(tag(ind))*xobs(:,ind));
        else
            col = cos(imag(tag(ind))*xobs(:,ind));
        end
    end
end