function A = seir_A(u,p,nP,it)
% -----------------------------------------------------------------------------
% diego domenzain
% spring 2020 while at Colorado School of Mines
% ------------------------------------------------------------------------------
% inversion of the SEIR model using the adjoint method.
% this function computes the A matrix.
% ------------------------------------------------------------------------------
b = p(1);                % infection rate
a = p(2);                % protection rate 
g = p(3);                % inverse of average latent time
d = p(4);                % inverse of average quarantine time
l = p(5);                % cure rate
k = p(6);                % mortality 
% ------------------------------------------------------------------------------
A      = zeros(7,7);
% ------------------------------------------------------------------------------
A(1,1) = -(b*u(it,3)/nP) - a;
A(2,1) =  b*u(it,3)/nP;
A(7,1) =  a;
A(2,2) = -g;
A(3,2) =  g;
A(3,3) = -d;
A(4,3) =  d;
A(4,4) = -(l+k);
A(5,4) =  l;
A(6,4) =  k;
% ------------------------------------------------------------------------------
end