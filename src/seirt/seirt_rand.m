function p = seirt_rand(corona_)

E = Inf;
for i_=1:corona_.nrand
 corona_.p = rand(8,1);
 corona_   = seirt_inv(corona_);
 if corona_.E<E
  E = corona_.E;
  p = corona_.p;
 end
end

end