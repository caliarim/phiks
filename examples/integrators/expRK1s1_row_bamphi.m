function [U,info] = expRK1s1_row_bamphi(tstar,ts,Afun,U,g)
% Method: exponential Euler

  kappa = 2^(-12); % adr
  t = 0;
  tau = tstar/ts;
  tol = kappa*tau^2;
  V(:,1) = U;
  info = [];
  opts.norm = 2;
  for j = 1:ts
    gtU = g(t,V(:,1));
    V(:,2) = gtU;
    normU = norm(V(:,1));
    opts.tol = tol*normU;
    [V(:,1),info] = bamphi(tau,Afun,[],V(:,1:2),opts,info);
    t = t+tau;
  end
  U = V(:,1);
