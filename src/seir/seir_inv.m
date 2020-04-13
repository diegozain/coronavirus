function corona_ = seir_inv(corona_)
% -----------------------------------------------------------------------------
% diego domenzain
% spring 2020 while at Colorado School of Mines
% ------------------------------------------------------------------------------
% inversion of the SEIR model using the adjoint method.
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
 corona_.u  = seir_fwd(corona_);
 % -----------------------------------------------------------------------------
 % error 
 corona_.d  = seir_u2d(corona_);
 corona_.err= corona_.d - corona_.do;
 corona_.E  = sum(sum(corona_.err.^2,2),1)*corona_.dt;
 % -----------------------------------------------------------------------------
 % adjoint and gradient
 corona_.adj  = seir_adj(corona_);
 corona_.grad_= seir_grad(corona_);
 % -----------------------------------------------------------------------------
 % regularization can go here
 % corona_.grad_ = corona_.grad_ .* [1e-2; 1e+1; 1e+1; 1e+1; 1e+1; 1e+1];
 % corona_.grad_ = corona_.grad_ + corona_.bo*(corona_.p-corona_.po);
 % -----------------------------------------------------------------------------
 corona_.update_ = - corona_.grad_;
 % ---------------------------------------------------------------------------
 % bfgs
 q        =  corona_.p - corona_.p_;
 y        =  corona_.grad_ - corona_.grad__;
 corona_.H=  broflegos(corona_.H,q,y);
 % ---------------------------------------------------------------------------
 corona_.grad__= corona_.grad_;
 corona_.p_    = corona_.p;
 % ---------------------------------------------------------------------------
 if ~isfield(corona_,'boH')
   boH = 1e-2*abs(mean(corona_.H(:)));
   % boH = max(boH,1e-6);
 else
   boH = corona_.boH;
 end
 corona_.update_  = -(corona_.H + boH*eye(size(corona_.H) )) \ corona_.grad_;
 % -----------------------------------------------------------------------------
 % variable step size
 if isfield(corona_,'step_')
   step_ = corona_.step_;
 elseif isfield(corona_,'perturb')
   step_ = seir_step(corona_);
 elseif isfield(corona_,'step_ave') 
   step_ = 1/5/norm(corona_.update_);
 else 
   step_ = 1;
 end
 % -----------------------------------------------------------------------------
 update_ = step_*corona_.update_;
 % -----------------------------------------------------------------------------
 % % momentum 
 % update_ = update_ + 2e-1*update__;
 % update__= update_;
 % -----------------------------------------------------------------------------
 % update
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