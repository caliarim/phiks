function [U,m_kry] = expRK1s1_row_phipmsimuliom(tstar,ts,Afun,U,g)
% Method: exponential Euler

  kappa = 2^(0); % adr
  t = 0;
  tau = tstar/ts;
  tol = kappa*tau^2;
  V(:,1) = U;
  m_kry = 1;
  iom = 2;
  for j = 1:ts
    gtU = g(t,V(:,1));
    V(:,2) = gtU;
    normU = norm(V(:,1));
    krytol = tol*normU;
    [V(:,1),m_kry] = phipm_simul_iom(tau,Afun,V(:,1:2),krytol,m_kry,iom);
    t = t+tau;
  end
  U = V(:,1);
