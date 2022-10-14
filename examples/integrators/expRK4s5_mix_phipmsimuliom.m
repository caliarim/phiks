function [U,m_kry] = expRK4s5_mix_phipmsimuliom2(tstar,ts,Afun,U,g)
% Method: exponential Runge Kutta of order 4

  if ts == 20
    kappa = 2^(-7); % brusselator
  else
    kappa = 2^(-10);
  end
  c(2)=1/2;
  c(3)=1/2;
  c(4)=1;
  c(5)=1/2;
  tau = tstar/ts;
  tol = kappa*tau^5;
  t = 0;
  m_kry1 = 1;
  m_kry2 = 1;
  m_kry3 = 1;
  m_kry4 = 1;
  m_kry5 = 1;
  iom = 2;

  for n = 1:ts
    gtU = g(t,U);
    normU = norm(U(:));
    krytol = tol*normU;
    V(:,1) = zeros(size(U));
    V(:,2) = Afun(U) + tau*gtU;
    [aux,m_kry1]=phipm_simul_iom([c(2),c(4)],Afun,V(:,1:2),krytol,m_kry1,iom);

    D2 = U + aux(:,1);
    D2 = g(t+c(2)*tau,D2)-gtU;

    V(:,2) = zeros(size(U));
    V(:,3) = tau*D2;
    [aux2,m_kry2]=phipm_simul_iom([c(3),c(4)],Afun,V(:,1:3),krytol,m_kry2,iom);

    D3 = U + aux(:,1) + aux2(:,1)/(c(3)^2);
    D3 = g(t+c(3)*tau,D3)-gtU;

    V(:,3) = tau*D3;
    [aux3,m_kry3]=phipm_simul_iom([c(5),c(4)],Afun,V(:,1:3),krytol,m_kry3,iom);

    D4 = U + aux(:,2) + aux2(:,2)/(c(4)^2) + aux3(:,2)/(c(4)^2);
    D4 = g(t+c(4)*tau,D4)-gtU;

    V(:,3) = -1/4*tau*D4;
    V(:,4) = tau*(-D2-D3+D4);
    [aux4,m_kry4] = phipm_simul_iom([c(5),c(4)],Afun,V,krytol,m_kry4,iom);

    D5 = U + aux(:,1) + 1/2*aux2(:,1)/(c(3)^2) + 1/4*aux2(:,2)/(c(4)^2) + ...
         1/2*aux3(:,1)/(c(5)^2) + 1/4*aux3(:,2)/(c(4)^2) +  ...
         4*aux4(:,1) + aux4(:,2);
    D5 = g(t+c(5)*tau,D5)-gtU;

    V(:,3) = tau*(-D4+4*D5);
    V(:,4) = tau*(4*D4-8*D5);
    [Utmp,m_kry5] = phipm_simul_iom(1,Afun,V,krytol,m_kry5,iom);
    U = U + aux(:,2) + Utmp;
    t = t+tau;
  end
  m_kry = [m_kry1,m_kry2,m_kry3,m_kry4,m_kry5];
