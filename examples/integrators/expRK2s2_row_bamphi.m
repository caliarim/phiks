function [U,info] = expRK2s2_row_bamphi(tstar,ts,Afun,U,g)
% Method: ETD2RK

  kappa = 2^(-9); % adr
  t = 0;
  tau = tstar/ts;
  tol = kappa*tau^3;
  V(:,1) = U;
  info1 = [];
  info2 = [];
  opt.norm = 2;
  for j = 1:ts
    gtU = g(t,V(:,1));
    V(:,2) = gtU;
    normU = norm(V(:,1));
    opts.tol = tol*normU;
    [U2,info1] = bamphi(tau,Afun,[],V(:,1:2),opts,info1);
    V(:,3) = g(t+tau,U2)-gtU;
    [V(:,1),info2] = bamphi(tau,Afun,[],V/diag([1,1,tau]),opts,info2);
    t = t+tau;
  end
  U = V(:,1);
  info = {info1,info2};
