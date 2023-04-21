% f = @(x,y)x+y;
% h = @(x,y)x.*y;
% fcns1 = {f,h};
% 
% g = @(x,y)x-y;
% q = @(x,y)x.^y; 
% fcns2 = {g,q};
% 
% nstates = 2;
% prodfcn_cell = prodlib(nstates,fcns1,fcns2);
% cellfun(@(f) f(2,3), prodfcn_cell)

function prodfcn_cell = prodlib(nstates,fcns1,fcns2)
    J1 = length(fcns1);
    J2 = length(fcns2);
    prodfcn_cell = cell(J1,J2);
    xstr=reshape(strcat('x',num2str((1:nstates)'),',')',[],1)';
    xstr = xstr(1:end-1);
    for i=1:J1
        f = fcns1{i};
        for j=1:J2        
            g = fcns2{j};
%             prodfcn_cell{i,j} = eval(['@(',xstr,')f(',xstr,').*g(',xstr,')']);
            a = strfind(functions(g).function,')');
            b = strfind(functions(g).function,')');
            prodfcn_cell{i,j} = eval([functions(f).function(1:b),'(',...
                functions(f).function(b+1:end),').*(',functions(g).function(a(1)+1:end),')']);
%             setfield(functions(prodfcn_cell{i,j}),'function',[functions(f).function,functions(g).function])
        end
    end
    prodfcn_cell = prodfcn_cell(:); 
end