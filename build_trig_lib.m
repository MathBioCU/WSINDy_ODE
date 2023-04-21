% nstates = 3;
% trigs = [0:5];
% custom_tags = [0 0 1i*6];
% [poly_cell,tags] = build_trig_lib(nstates,polys,custom_tags);

function [trig_cell,tags] = build_trig_lib(nstates,trigs,custom_tags,filters)
    if and(~isempty(custom_tags),size(custom_tags,2)~=nstates)
        disp(['invalid custom poly tags'])
        return
    end
    
    tags = custom_tags;
    for k=1:length(trigs)
        trig_inds = [-trigs(k)*1i*eye(nstates);trigs(k)*1i*eye(nstates)];
        tags = [tags; trig_inds];
    end
    tags = unique(tags,'rows');
    for j=1:length(filters)
        inds = filters{j}(tags);
        tags = tags(inds,:);
    end
    J = size(tags,1);

    trig_cell = cell(J,1);
    xstr=reshape(strcat('x',num2str((1:nstates)'),',')',[],1)';
    xstr = xstr(1:end-1);
    for j=1:J
        ind = find(tags(j,:));
        if isempty(ind)
            k=0;
            ind=1;
        else
            k=imag(tags(j,ind));
        end
        if k<0            
            trig_cell{j} = eval(['@(',xstr,')sin(',num2str(-k),'*x',num2str(ind),')']);
        else
            trig_cell{j} = eval(['@(',xstr,')cos(',num2str(k),'*x',num2str(ind),')']);
        end
    end
end