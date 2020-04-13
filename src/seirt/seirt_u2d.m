function d = seirt_u2d(corona_)
d = zeros(corona_.nt,4);
for it=1:corona_.nt
 d_ = corona_.Mu2d * (corona_.u(it,:)).' ;
 d(it,:) = d_.';
end
end