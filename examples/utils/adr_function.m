function [A,Kfun,g,gv,U0,U0v,h,K] = adr_function(epsilon,alpha,n)
%
% function [A,Kfun,g,gv,U0,U0v,h,K] = adr_function(epsilon,alpha,n)
%
% A is the cell of (full) uni-directional matrices,
% Kfun is the function which performs the action of their Kronecker sum,
% g is the nonlinear function in tensor form,
% gv is the nonlinear function in vector form,
% U0 the initial data in tensor form,
% gv is the initial data in vector form,
% h the vector of spacial step-sizes,
% K is the Kronecker sum of the matrices in A.

  d = length(n);
  for mu = 1:d
    x{mu} = linspace(0,1,n(mu)+2)';
    x{mu} = x{mu}(2:n(mu)+1);
    h(mu) = 1/(n(mu)+1);
    D2{mu} = spdiags(ones(n(mu),1)*[1,-2,1]/h(mu)^2,-1:1,n(mu),n(mu));
    D1{mu} = spdiags(ones(n(mu),1)*[-1,1]/2/h(mu),[-1,1],n(mu),n(mu));
    A{mu} = epsilon*D2{mu}+alpha(mu)*D1{mu};
  end
  if (d == 1)
    K = A{1};
  else
    K = kronsum(A);
  end
  Kfun = @(u) K*u;
  for mu = 1:d
    A{mu} = full(A{mu});
  end
  [X{1:d}] = ndgrid(x{1:d});
  U0 = X{1}.*(1-X{1})*4;
  for mu = 2:d
    U0 = U0 .* X{mu}.*(1-X{mu})*4;
  end
  U0v = U0(:);
  Phi1 = zeros(size(X{1}));
  Phi2 = zeros(size(X{1}));
  for mu = 1:d
    aux1 = alpha(mu)*(1-2*X{mu})*4;
    aux2 = epsilon*(-2)*4;
    for eta = [1:mu-1,mu+1:d]
      aux1 = aux1.*X{eta}.*(1-X{eta})*4;
      aux2 = aux2.*X{eta}.*(1-X{eta})*4;
    end
    Phi1 = Phi1+aux1;
    Phi2 = Phi2+aux2;
  end
  g = @(t,U) 1./(1+U.^2)+(U0-Phi2-Phi1)*exp(t)-1./(1+U0.^2*exp(2*t));
  gv = @(t,U) 1./(1+U(:).^2)+(U0(:)-Phi2(:)-Phi1(:))*exp(t)-...
       1./(1+U0(:).^2*exp(2*t));
