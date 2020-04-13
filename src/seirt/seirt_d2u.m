function u = seirt_d2u(corona_)
u = zeros(corona_.nt,7);
for it=1:corona_.nt
 u(it,:) = corona_.d(it,:) * corona_.M.';
end
end