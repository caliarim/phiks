function [U,m_kry] = expRK2s2_row_phipmsimuliom(tstar,ts,Afun,U,g)
% Method: ETD2RK

  kappa = 2^(2); % adr
  t = 0;
  tau = tstar/ts;
  tol = kappa*tau^3;
  V(:,1) = U;
  m_kry1 = 1;
  m_kry2 = 1;
  iom = 2;
  for j = 1:ts
    gtU = g(t,V(:,1));
    V(:,2) = gtU;
    normU = norm(V(:,1));
    krytol = tol*normU;
    [U2,m_kry1] = phipm_simul_iom(tau,Afun,V(:,1:2),krytol,m_kry1,iom);
    V(:,3) = g(t+tau,U2)-gtU;
    [V(:,1),m_kry2] = phipm_simul_iom(tau,Afun,V/diag([1,1,tau]),krytol,m_kry2,iom);
    t = t+tau;
  end
  U = V(:,1);
  m_kry = [m_kry1,m_kry2];
