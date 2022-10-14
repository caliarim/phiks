function [U,info] = expRK3s3_column_bamphi(tstar,ts,Afun,U,g);
% Method: exponential Runge Kutta of order 3

  kappa = 2^4; % allencahn
  c(2)=1/4;
  c(3)=1/2;
  gamma = (3*c(3)-2)*c(3)/(2-3*c(2))/c(2);
  tau = tstar / ts;
  tol = kappa*tau^4;
  t = 0;
  info1 = [];
  info2 = [];
  info3 = [];
  opts.norm = 2;
  for j = 1:ts
    normU = norm(U(:));
    opts.tol = tol*normU;
    gtU = g(t,U);
    FtU(:,2) = Afun(U)+gtU;
    [aux1,info1] = bamphi(tau*[c(2),c(3),1],Afun,[],FtU,opts,info1);
    U2 = U+aux1(:,1);
    D2(:,3) = g(t+c(2)*tau,U2)-gtU;
    [aux2,info2] = bamphi(tau*[c(2),c(3),1],Afun,[],D2/tau,opts,info2);
    U3 = U+aux1(:,2)+(gamma/c(2))*aux2(:,1)+(1/c(2))*aux2(:,2);
    D3(:,3) = g(t+c(3)*tau,U3)-gtU;
    [aux3,info3] = bamphi(tau,Afun,[],D3/tau,opts,info3);
    U = U+aux1(:,3)+gamma/(gamma*c(2)+c(3))*aux2(:,3)+1/(gamma*c(2)+c(3))*aux3;
    t = t+tau;
  end
  info = {info1,info2,info3};
