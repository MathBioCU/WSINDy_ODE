function out=arrayfunvec(arr_in,f,d)
    dims=size(arr_in);
    dim = length(dims);
    e = ones(dim);
    e(1:dim+1:end)=dims;
    e=num2cell(e);
    inds=arrayfun(@(x)ones(e{x:dim:end}),1:dim,'uni',0);
    inds{d} = dims(d);
    arr_in_cell = mat2cell(arr_in,inds{:});
    out = cell2mat(cellfun(@(x)f(x),arr_in_cell,'uni',0));
end