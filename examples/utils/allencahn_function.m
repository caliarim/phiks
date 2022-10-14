function [A,Kfun,g,gv,U0,U0v,h] = allencahn_function(epsilon,N,alpha,n,FW)
%
% function [A,Kfun,g,gv,U0,U0v,h] = allencahn_function(epsilon,alpha,n,FW)
%
% A is the cell of (full) uni-directional matrices,
% Kfun is the function which performs the action of their Kronecker sum,
% g is the nonlinear function in tensor form,
% gv is the nonlinear function in vector form,
% U0 the initial data in tensor form,
% gv is the initial data in vector form,
% h the vector of spacial step-sizes.

  include_lin = 1;
  d = length(n);
  for mu = 1:d
    x{mu} = linspace(0,1,n(mu))';
    h(mu) = 1/(n(mu)-1);
    switch FW
      case 1 % FD
        A{mu} = spdiags(ones(n(mu),1)*[1,-2,1]/h(mu)^2,-1:1,n(mu),n(mu));
        A{mu}(1,1:2) = [-2,2]/h(mu)^2;
        A{mu}(n(mu),n(mu)-1:n(mu)) = [2,-2]/h(mu)^2;
        A{mu} = A{mu}+include_lin/(d*epsilon^2)*speye(n(mu));
      case 2 % FEM with lump (which coincides with FD)
        M{mu} = spdiags([h(mu)/2;h(mu)*ones(n(mu)-2,1);h(mu)/2],0,n(mu),n(mu));
        S{mu} = spdiags(ones(n(mu),1)*[1,-2,1]/h(mu),-1:1,n(mu),n(mu));
        S{mu}(1,1:2) = [-1,1]/h(mu);
        S{mu}(n(mu),n(mu)-1:n(mu)) = [1,-1]/h(mu);
        A{mu} = M{mu}\S{mu}+include_lin/(d*epsilon^2)*speye(n(mu));
      case 3 % FEM exact
        M{mu} = spdiags(ones(n(mu),1)*[1/6,2/3,1/6]*h(mu),-1:1,n(mu),n(mu));
        M{mu}(1,1:2) = [1/3,1/6]*h(mu);
        M{mu}(n(mu),n(mu)-1:n(mu)) = [1/6,1/3]*h(mu);
        S{mu} = spdiags(ones(n(mu),1)*[1,-2,1]/h(mu),-1:1,n(mu),n(mu));
        S{mu}(1,1:2) = [-1,1]/h(mu);
        S{mu}(n(mu),n(mu)-1:n(mu)) = [1,-1]/h(mu);
        A{mu} = M{mu}\S{mu}+include_lin/(d*epsilon^2)*speye(n(mu));
      case 4 % FD4 with hom. Neumann
        A{mu} = spdiags(ones(n(mu),1)*[-1/12,4/3,-5/2,4/3,-1/12]/h(mu)^2,...
                           -2:2,n(mu),n(mu));
        A{mu}(1,1:5) = [-145/3,56,-6,-8/3,1]/(12*h(mu)^2);
        A{mu}(2,1:4) = [29/18,-3,3/2,-1/9]/h(mu)^2;
        A{mu}(n(mu),n(mu)-4:n(mu)) = [1,-8/3,-6,56,-145/3]/(12*h(mu)^2);
        A{mu}(n(mu)-1,n(mu)-3:n(mu)) = [-1/9,3/2,-3,29/18]/h(mu)^2;
        A{mu} = A{mu}+include_lin/(d*epsilon^2)*speye(n(mu));
      case 5 % FD + different way to treat linearity
        M{mu} = spdiags([h(mu)/2;h(mu)*ones(n(mu)-2,1);h(mu)/2],0,n(mu),n(mu));
        S{mu} = spdiags(ones(n(mu),1)*[1,-2,1]/h(mu),-1:1,n(mu),n(mu));
        S{mu}(1,1:2) = [-1,1]/h(mu);
        S{mu}(n(mu),n(mu)-1:n(mu)) = [1,-1]/h(mu);
        A{mu} = M{mu}\S{mu}+epsilon^2/d*speye(n(mu));
      otherwise
        error('FW can be just 1,2,5 (FD2/lumped FEM), 3 (FEM) or 4 (FD4)')
    end
  end
  if (d == 1)
    K = A{1};
    Kfun = @(u) K*u;
  else
    switch FW
      case {1,2,4,5} % faster sparse matvec product
        K = kronsum(A);
        Kfun = @(u) K*u;
        for mu = 1:d
          A{mu} = full(A{mu});
        end
      case 3
        for mu = 1:d
          A{mu} = full(A{mu});
        end
        Kfun = @(u) reshape(kronsumv(reshape(u,n),A),[],1);
    end
  end
  [X{1:d}] = ndgrid(x{1:d});
  %theta = atan((X{2}-0.5)./(X{1}-0.5));
  %theta = pi*(X{1}<=0.5)+theta;
  theta = atan2(X{2}-0.5,X{1}-0.5);
  U0 = tanh((0.25+0.1*cos(N*theta)-sqrt((X{1}-0.5).^2+(X{2}-0.5).^2))/...
            (alpha*sqrt(2)));
  U0v = U0(:);
  if FW == 5
    g = @(t,U) 1/epsilon^2*(U-U.^3)-epsilon^2*U;
    gv = @(t,U) 1/epsilon^2*(U(:)-U(:).^3)-epsilon^2*U(:);
  else
    g = @(t,U) (1-include_lin)/epsilon^2*U-1/epsilon^2*U.^3;
    gv = @(t,U) -1/epsilon^2*U(:).^3;
  end
