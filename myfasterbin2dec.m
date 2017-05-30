function y = myfasterbin2dec(table)

%
% table: should be a vector (i.e. of size Nx1 or 1xN)
%

[m n]   = size(table);
if min(m,n) > 1
    error('please input a vector in myfasterbin2dec(.)');
end
if max(table)>1
    error('input is not binary');
end
if m > n
    table = table';
end
L       = length(table);
expo    = (L-1):-1:0;
vect    = (2.^expo)';
y       = table*vect;

end