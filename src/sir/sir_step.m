function step_ = sir_step(corona_)
% -----------------------------------------------------------------------------
% diego domenzain
% spring 2020 while at Colorado School of Mines
% ------------------------------------------------------------------------------
% inversion of the SEIR model using the adjoint method.
% this function performs a parabolic search for the optimum step size
% ------------------------------------------------------------------------------
E       = corona_.E;
perturb = corona_.perturb;
grad_   = corona_.grad_;
% grad_   = normali(grad_);
p       = corona_.p;
% ------------------------------------------------------------------------------
E_ = [];
for i_=1:numel(perturb)
  % ----------------------------------------------------------------------------
  % perturb 
  corona_.p  = p.*exp(p.*(-perturb(i_)*grad_));
  % ----------------------------------------------------------------------------
  % forward model
  corona_.f  = zeros(numel(corona_.t),1);
  u          = sir_fwd(corona_);
  err        = u - corona_.uo;
  E__        = sum(sum(err.^2,2),1)*corona_.dt;
  % ----------------------------------------------------------------------------
  E_ = [E_ ; E__];
end
% ------------------------------------------------------------------------------
% bundle together
E = [E ; E_];
k = [0 ; perturb];
% ------------------------------------------------------------------------------
% parabola approx
warning off;
k = polyfit(k,E,2);
warning on;
% find zero of parabola (update = -step*gradient)
step_ = -k(2)/(2*k(1));
% ------------------------------------------------------------------------------
end