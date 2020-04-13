function grad_ = seirt_grad(corona_)
% ------------------------------------------------------------------------------
% diego domenzain
% spring 2020 while at Colorado School of Mines
% ------------------------------------------------------------------------------
% inversion of the SEIR model using the adjoint method.
% this function computes the gradient of the objective function SSR.
% ------------------------------------------------------------------------------
u  = corona_.u;
p  = corona_.p;
t  = corona_.t;
nP = corona_.nP;
adj= corona_.adj;
dt = corona_.dt;
nt = numel(corona_.t);
% ------------------------------------------------------------------------------
grad_= zeros(1,numel(corona_.p));
% ------------------------------------------------------------------------------
for it=1:nt
 it_  = nt-(it-1);
 % -----------------------------------------------------------------------------
 Apu  =  seirt_Apu(u,p,t,nP,it_);
 % -----------------------------------------------------------------------------
 adj_ = adj(it,:);
 % -----------------------------------------------------------------------------
 grad_= grad_ + (adj_*Apu); 
 % -----------------------------------------------------------------------------
 % regularization goes here
 % grad_ = grad_ + regularization term ;
 % -----------------------------------------------------------------------------
end
grad_ = grad_*dt;
grad_ = grad_.';
end