function [A,Kfun,g,gv,U0,U0v,h,M,S] = brusselator_function(a,b,c,d1,d2,n,FW)
%
% function [A,Kfun,g,gv,U0,U0v,h,M,S] = brusselator_function(a,b,c,d1,d2,n,FW)
%
% A is the cell of (full) uni-directional matrices,
% Kfun is the function which performs the action of their Kronecker sum,
% g is the nonlinear function in tensor form,
% gv is the nonlinear function in vector form,
% U0 the initial data in tensor form,
% gv is the initial data in vector form,
% h the vector of spacial step-sizes,
% M are the uni-directional mass matrices (if FEM is chosen),
% S are the uni-directional stiffness matrices (if FEM is chosen).

  d = length(n);
  prodn = prod(n);
  for mu = 1:d
    x{mu} = linspace(0,1,n(mu))';
    h(mu) = 1/(n(mu)-1);
    switch FW
      case {1,2} % FD or FEM with lump (they coincide) with hom. Neumann
        A{1}{mu} = spdiags(ones(n(mu),1)*[1,-2,1]/h(mu)^2,-1:1,n(mu),n(mu));
        A{1}{mu}(1,1:2) = [-2,2]/h(mu)^2;
        A{1}{mu}(n(mu),n(mu)-1:n(mu)) = [2,-2]/h(mu)^2;
        A{2}{mu} = d2*A{1}{mu};
        A{1}{mu} = d1*A{1}{mu}-(b+1)/d*speye(n(mu));
      case 3 % FEM exact with hom. Neumann
        M{mu} = spdiags(ones(n(mu),1)*[1/6,2/3,1/6]*h(mu),-1:1,n(mu),n(mu));
        M{mu}(1,1:2) = [1/3,1/6]*h(mu);
        M{mu}(n(mu),n(mu)-1:n(mu)) = [1/6,1/3]*h(mu);
        S{mu} = spdiags(ones(n(mu),1)*[1,-2,1]/h(mu),-1:1,n(mu),n(mu));
        S{mu}(1,1:2) = [-1,1]/h(mu);
        S{mu}(n(mu),n(mu)-1:n(mu)) = [1,-1]/h(mu);
        A{2}{mu} = d2*(M{mu}\S{mu});
        A{1}{mu} = d1/d2*(A{2}{mu})-(b+1)/d*speye(n(mu));
      case 4 % FD4 with hom. Neumann
        A{1}{mu} = spdiags(ones(n(mu),1)*[-1/12,4/3,-5/2,4/3,-1/12]/h(mu)^2,...
                           -2:2,n(mu),n(mu));
        A{1}{mu}(1,1:5) = [-145/3,56,-6,-8/3,1]/(12*h(mu)^2);
        A{1}{mu}(2,1:4) = [29/18,-3,3/2,-1/9]/h(mu)^2;
        A{1}{mu}(n(mu),n(mu)-4:n(mu)) = [1,-8/3,-6,56,-145/3]/(12*h(mu)^2);
        A{1}{mu}(n(mu)-1,n(mu)-3:n(mu)) = [-1/9,3/2,-3,29/18]/h(mu)^2;
        A{2}{mu} = d2*A{1}{mu};
        A{1}{mu} = d1*A{1}{mu}-(b+1)/d*speye(n(mu));
      otherwise
        error('FW can take as input 1,2 (FD/lumped FEM) or 3 (FEM) or 4 (FD4)')
    end
  end
  if (d == 1)
    Kfun = @(u) [A{1}{1}*u(1:n(1));A{2}{1}*u(n(1)+1:2*n(1))];
  else
    switch FW
      case {1,2,4} % faster sparse matvec product
        K{1} = kronsum(A{1});
        K{2} = kronsum(A{2});
        Kfun = @(u) [K{1}*u(1:prodn);...
                     K{2}*u(prodn+1:2*prodn)];
        for mu = 1:d
          A{1}{mu} = full(A{1}{mu});
          A{2}{mu} = full(A{2}{mu});
        end
      case 3
        for mu = 1:d
          A{1}{mu} = full(A{1}{mu});
          A{2}{mu} = full(A{2}{mu});
        end
        Kfun = @(u)[reshape(kronsumv(reshape(u(1:prodn),n),A{1}),[],1);...
                   reshape(kronsumv(reshape(u(prodn+1:2*prodn),n),A{2}),[],1)];
    end
  end
  [X{1:d}] = ndgrid(x{1:d});
  U0{1} = 16 * (X{1}.*(1-X{1})).^2;
  for mu = 2:d
    U0{1} = 16 * U0{1}.*(X{mu}.*(1-X{mu})).^2;
  end
  U0{2} = c*ones(n);
  U0v = [U0{1}(:);U0{2}(:)];
  g{1} = @(t,U,V) a+U.^2.*V;
  g{2} = @(t,U,V) -U.^2.*V+b*U;
  gv = @(t,u) [a+u(1:prodn).^2.*u(prodn+1:2*prodn);...
               -u(1:prodn).^2.*u(prodn+1:2*prodn)+b*u(1:prodn)];
