function sumfcn_cell = sumlib(nstates,fcns1,fcns2)
    J1 = length(fcns1);
    J2 = length(fcns2);
    sumfcn_cell = cell(J1,J2);
    xstr=reshape(strcat('x',num2str((1:nstates)'),',')',[],1)';
    xstr = xstr(1:end-1);
    for i=1:J1
        f = fcns1{i};
        for j=1:J2        
            g = fcns2{j};
            sumfcn_cell{i,j} = eval(['@(',xstr,')f(',xstr,')+g(',xstr,')']);
        end
    end
    sumfcn_cell = sumfcn_cell(:); 
end