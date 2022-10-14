function U = expRK3s3_column_sylvphi(tstar,ts,A,U,g)
% Method: exponential Runge Kutta of order 3

  c(2) = 1/4;
  c(3) = 1/2;
  gamma = ((3*c(3)-2)*c(3))/((2-3*c(2))*c(2));
  tau = tstar/ts;
  t = 0;

  Ac2{1} = c(2)*tau*A{1};
  Ac2{2} = c(2)*tau*A{2};
  Ac3{1} = 2*Ac2{1};
  Ac3{2} = 2*Ac2{2};
  A1{1} = 2*Ac3{1};
  A1{2} = 2*Ac3{2};

  for mu = 1:2
    Ac2Phi1Ac2{mu} = Ac2{mu}*phipade(Ac2{mu},1);
    Ac3Phi1Ac3{mu} = Ac3{mu}*phipade(Ac3{mu},1);
    Ac22Phi2Ac2{mu} = Ac2{mu}*(Ac2{mu}*phipade(Ac2{mu},2));
    Ac32Phi2Ac3{mu} = Ac3{mu}*(Ac3{mu}*phipade(Ac3{mu},2));
    A1Phi1A1{mu} = A1{mu}*phipade(A1{mu},1);
    A12Phi2A1{mu} = A1{mu}*(A1{mu}*phipade(A1{mu},2));
  end

  for j = 1:ts
    F = A{1}*U+U*A{2}.' + g(t,U);

    U2 = sylvphi(Ac2,F,1,Ac2Phi1Ac2);
    U2 = U + c(2)*tau*U2;
    D2 = g(t+c(2)*tau,U2) - g(t,U);

    U3 = sylvphi(Ac3,F,1,Ac3Phi1Ac3); %phi1(1/2)F

    U3a = sylvphi(Ac2,D2,2,Ac2Phi1Ac2,Ac22Phi2Ac2); %phi2(1/4)D2

    U3b = sylvphi(Ac3,D2,2,Ac3Phi1Ac3,Ac32Phi2Ac3); %phi2(1/2)D2

    U3 = U + c(3)*tau*U3 + gamma*c(2)*tau*U3a + c(3)^2/c(2)*tau*U3b;
    D3 = g(t+c(3)*tau,U3) - g(t,U);

    Ua = sylvphi(A1,F,1,A1Phi1A1); %phi1(1)F
    D3 = gamma/(gamma*c(2)+c(3))*D2 + 1/(gamma*c(2)+c(3))*D3;

    Ub = sylvphi(A1,D3,2,A1Phi1A1,A12Phi2A1); %phi2(1)D3

    U = U + tau*Ua + tau*Ub;
    t = t + tau;
  end
