function corona_ = seirt_inv(corona_)
% -----------------------------------------------------------------------------
% diego domenzain
% spring 2020 while at Colorado School of Mines
% ------------------------------------------------------------------------------
% inversion of the SEIRt model using the adjoint method.
% 
% ------------------------------------------------------------------------------
iter = 0;
corona_.E = Inf;
update__  = zeros(size(corona_.p));
E_history = [];
s_history = [];
p_history = corona_.p;
fprintf('\n\n               the inversion has begun! \n\n');
% ------------------------------------------------------------------------------
while (iter<corona_.niter & corona_.E>corona_.tol_error)
 % -----------------------------------------------------------------------------
 % forward model
 corona_.f  = zeros(numel(corona_.t),1);
 corona_.u  = seirt_fwd(corona_);
 % -----------------------------------------------------------------------------
 % error 
 corona_.d  = seirt_u2d(corona_);
 corona_.err= corona_.d - corona_.do;
 corona_.E  = sum(sum(corona_.err.^2,2),1)*corona_.dt;
 % -----------------------------------------------------------------------------
 % adjoint and gradient
 corona_.adj  = seirt_adj(corona_);
 corona_.grad_= seirt_grad(corona_);
 % -----------------------------------------------------------------------------
 % regularization can go here
 
 % -----------------------------------------------------------------------------
 % store the gradient as update so the step-size computation is independent of 
 % grad descent or gauss-newton routines
 corona_.update_ = - corona_.grad_;
 % ---------------------------------------------------------------------------
 % 
 %                        second order methods
 % 
 % ---------------------------------------------------------------------------
 q =  corona_.p - corona_.p_;
 y =  corona_.grad_ - corona_.grad__;
 % ---------------------------------------------------------------------------
 %                            bfgs = Hessian
 % ---------------------------------------------------------------------------
 % approximate the hessian
 corona_.H        =  broflegos(corona_.H,q,y);
 % ---------------------------------------------------------------------------
 % find boh
 
 % ---------------------------------------------------------------------------
 % bgfs update
 corona_.update_  = -(corona_.H + boH*eye(size(corona_.H) )) \ corona_.grad_;
 % ---------------------------------------------------------------------------
 %                            bfgs = inverse( Hessian )
 % ---------------------------------------------------------------------------
 % approximate the inverse of the hessian
 corona_.B      =  broflegos(corona_.B,q,y);
 % ---------------------------------------------------------------------------
 % bgfs update
 corona_.update_= - corona_.B * corona_.grad_;
 % ---------------------------------------------------------------------------
 % store gradient and parameters for next iteration
 corona_.grad__= corona_.grad_;
 corona_.p_    = corona_.p;
 % -----------------------------------------------------------------------------
 % 
 %                     find the right step-size
 % 
 % -----------------------------------------------------------------------------
 step_ = seirt_step(corona_);
 % -----------------------------------------------------------------------------
 % compute the update with the right step-size 
 update_ = step_*corona_.update_;
 % -----------------------------------------------------------------------------
 % momentum 
 
 % -----------------------------------------------------------------------------
 %                                 update
 % -----------------------------------------------------------------------------
 corona_.p  = corona_.p.*exp(corona_.p.*(update_));
 % -----------------------------------------------------------------------------
 E_history  = [E_history corona_.E];
 % -----------------------------------------------------------------------------
 iter = iter+1;
 % print some funny stuff so user is amused
 if mod(iter,fix(corona_.niter*0.1))==1
   fprintf('  current error %2.2d at iteration %i \n',corona_.E,iter);
   p_history  = [p_history abs(update_)];
   s_history  = [s_history step_];
 end
end
% ------------------------------------------------------------------------------
fprintf('\n       total iterations were %i \n',iter);
corona_.p_history = p_history;
corona_.s_history = s_history;
corona_.E_history = E_history;
% ------------------------------------------------------------------------------
end