% nstates = 3;
% polys = [0 1 3];
% custom_tags = [0 0 4];
% [poly_cell,tags] = build_poly_lib(nstates,polys,custom_tags);

function [poly_cell,tags] = build_poly_lib(nstates,polys,custom_tags,filters)
    if and(~isempty(custom_tags),size(custom_tags,2)~=nstates)
        disp(['invalid custom poly tags'])
        return
    end
    
    tags = custom_tags;
    for p = 1:length(polys)
        monom_powers = partitionNk(polys(p),nstates);
        tags = [tags;monom_powers];
    end
    tags = unique(tags,'rows');
    for j=1:length(filters)
        inds = filters{j}(tags);
        tags = tags(inds,:);
    end
    J = size(tags,1);

    poly_cell = cell(J,1);
   
    xstr=reshape(strcat('x',num2str((1:nstates)'),',')',[],1)';
    xstr = xstr(1:end-1);
    for j=1:J
        xstr2=arrayfun(@(x,y) ['x',x,'.^(',y,').*'],num2str((1:nstates)'),num2str(tags(j,:)'),'uni',0);
        xstr2 = strcat(xstr2{:});
        xstr2 = xstr2(1:end-2);
        poly_cell{j} = eval(['@(',xstr,')',xstr2]);
    end
end
% 
% initial symbolic version:
% function [poly_cell,tags] = build_poly_lib(nstates,polys,custom_tags)
%     if and(~isempty(custom_tags),size(custom_tags,2)~=nstates)
%         disp(['invalid custom poly tags'])
%         return
%     end
%     
%     tags = custom_tags;
%     for p = 1:length(polys)
%         monom_powers = partitionNk(polys(p),nstates);
%         tags = [tags;monom_powers];
%     end
%     tags = unique(tags,'rows');
%     J = size(tags,1);
% 
%     poly_cell = cell(J,1);
%    
%     X = str2sym(strcat('x',num2str((1:nstates)')));
%     Xc = sym2cell(X);
%     for j=1:J
%         if all(~tags(j,:))
%             xstr=reshape(strcat('x',num2str((1:nstates)'),',')',[],1)';
%             xstr = xstr(1:end-1);
%             poly_cell{j} = eval(['@(',xstr,')x1*0+1']);
%         else
%             poly_cell{j} = matlabFunction(prod(X.^(tags(j,:)')),'vars',Xc);
%         end
%     end
% end