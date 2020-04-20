function it = find_zero(u,val)
% find index of first 'val' crossing in vector u
it=1;
while u(it)>val
 it=it+1;
end
end