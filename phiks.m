function [varargout] = phiks(tau,A,V,p,tol,shat,shift)
%PHIKS phi-functions of a Kronecker sum applied to tensors.
%   PHIKS(TAU,A,V,P,TOL) when V is a cell of P+1 ND-arrays of the same size,
%   returns the linear combination
%
%   EXP(TAU*K)*V{1}+PHI_1(TAU*K)*V{2}+...+PHI_P(TAU*K)*V{P+1},
%
%   where K is the Kronecker sum of A{D},A{D-1},...,A{1}, with A{MU} a full
%   square matrix of size N(MU)xN(MU). No matrix K is assembled.
%   It is possible to use the scalar 0 instead of a ND-array of zero
%   elements in all the entries of the cell V, but the last.
%   The output is a ND-array of the same size of V{1}.
%
%   [COMBPHI_1,COMBPHI_2,...,COMBPHI_SCALES] = PHIKS(TAU,A,V,P,TOL,SCALES)
%   returns the linear combinations
%
%   EXP(TAU/2^(J-1)*K)*V{1}+1/2^(J-1)*PHI_1(TAU/2^(J-1)*K)*V{2}+...
%                      +(1/2^(J-1))^P*PHI_P(TAU/2^(J-1)*K)*V{P+1}.
%
%   for J=1,2,...,SCALES.
%
%   PHIKS(TAU,A,V,P,TOL) when V is a single ND-array, returns a cell of
%   ND-arrays corresponding to
%
%   PHI_{ELL-1}(TAU*K)*V,  ELL=1,2,...,P+1.
%
%   [PHI_1,PHI_2,...,PHI_SCALES] = PHIKS(TAU,A,V,P,TOL,SCALES) returns
%
%   PHI_{ELL-1}(TAU/2^(J-1)*K)*V, ELL=1,2,...,P+1.
%
%   for J=1,2,...,SCALES.
%
%   [___] = PHIKS(___, FALSE) disables the trace shifting strategy.
%
%   [___, S, Q, C] = PHIKS(___)  also returns the used scaling S, the number
%   employed of quadrature points Q, and the cost C in terms of Tucker
%   operators or equivalent operations.
  d = length(A);
  if (nargin <= 6)
    shift = true;
  end
  if (nargin <= 5)
    shat = 1;
  end
  if (p == 0) % require only exponential
    for mu = 1:d
      E{mu} = expm(tau/shat*A{mu});
    end
    if (iscell(V)) % linear combination phi
      varargout{shat} = tucker(V{1},E);
      for j = shat-1:-1:1
        varargout{j} = tucker(varargout{j+1},E);
      end
    else % phi same vector
      varargout{shat}{1} = tucker(V,E);
      for j = shat-1:-1:1
        varargout{j}{1} = tucker(varargout{j+1}{1},E);
      end
    end
    varargout{shat+1} = shat;
    varargout{shat+2} = 0;
    varargout{shat+3} = shat;
    return
  end
  % shift and rectangle
  for mu = 1:d
    de = diag(A{mu});
    na = length(de);
    sigma(mu) = shift*sum(de)/na;
    eigenA{mu} = rect(A{mu},'norm2');
    eigenA{mu} = tau*eigenA{mu};
    A{mu}(1:na+1:na*na) = de-sigma(mu);
  end
  eigentauK = sum(cat(1,eigenA{:}),1);
  Xi_boundary = fov_contour(eigentauK);
  if (iscell(V)) % linear combination phi
    zeroV = false (1, p+1);
    normV = zeros (1, p+1);
    if (isscalar(V{1}) && V{1} == 0) % find V{1} in fact the scalar 0
      zeroV(1) = true;
    end
    for ell = 1:p % find other V{ell} in fact the scalar 0
      if (isscalar(V{ell+1}) && V{ell+1} == 0)
        zeroV(ell+1) = true;
      else
        normV(ell+1) = norm(V{ell+1}(:),2);
      end
    end
    [s,q] = findsq_kernel_lcp(Xi_boundary,p,normV(2:p+1),tol,d,shat);
    [Psi,E,T] = integralgll_lcp(tau,A,V,p,s,q,sigma,zeroV);
    for j = s:-1:shat % squaring
      for ell = p:-1:1+(p-1)*(j == 1) % skip lower order combinationsw
        Psi{ell} = Psi{ell}+tucker(Psi{ell},E);
        factlminusk = 2^j;
        for k = ell-1:-1:1
          Psi{ell} = Psi{ell}+Psi{k}/factlminusk; % /2^((l-k)*j);
          factlminusk = factlminusk*(ell-k+1)*2^j;
        end
      end
      T = T+1+(p-1)*(j~=1);
      for mu = 1:d
        E{mu} = E{mu} * E{mu};
      end
    end
    for j = shat-1:-1:1 % squaring
      if (~zeroV(1))
        varargout{j+1} = tucker(V{1},E)+Psi{p};
        T = T+1;
      else
        varargout{j+1} = Psi{p};
      end
      for ell = p:-1:1+(p-1)*(j == 1) % skip lower order combinations
        Psi{ell} = Psi{ell}+tucker(Psi{ell},E);
        factlminusk = 2^j;
        for k = ell-1:-1:1
          Psi{ell} = Psi{ell}+Psi{k}/factlminusk; % /2^((l-k)*j);
          factlminusk = factlminusk*(ell-k+1)*2^j;
        end

      end
      T = T+1+(p-1)*(j~=1);
      for mu = 1:d
        E{mu} = E{mu} * E{mu};
      end
    end
    if (~zeroV(1))
      varargout{1} = tucker(V{1},E)+Psi{p};
      T = T+1;
    else
      varargout{1} = Psi{p};
    end
  else % phi same vector
    [s,q] = findsq_kernel_psv(Xi_boundary,p,norm(V(:),2),tol,d,shat);
    [varargout{shat},E,T] = integralgll_psv(tau,A,V,p,s,q,sigma);
    for j = s:-1:shat % squaring
      for ell = p:-1:1
        varargout{shat}{ell+1} = (varargout{shat}{ell+1}+...
                                    tucker(varargout{shat}{ell+1},E))/2^ell;
        factlminusk = 2^ell;
        for k = ell-1:-1:1
          varargout{shat}{ell+1} = varargout{shat}{ell+1}+...
                                     varargout{shat}{k+1}/factlminusk;
          factlminusk = factlminusk*(ell-k+1);
        end
      end
      for mu = 1:d
        E{mu} = E{mu} * E{mu};
      end
    end
    T = T+(s-shat+1)*p;
    if (s >= shat)
      varargout{shat}{1} = tucker (V, E);
      T = T+1;
    end
    for j = shat-1:-1:1
      for ell = p:-1:1
        varargout{j}{ell+1} = (varargout{j+1}{ell+1}+...
                               tucker(varargout{j+1}{ell+1},E))/2^ell;
        factlminusk = 2^ell;
        for k = ell-1:-1:1
          varargout{j}{ell+1} = varargout{j}{ell+1}+...
                                varargout{j+1}{k+1}/factlminusk;
          factlminusk = factlminusk*(ell-k+1);
        end
      end
      for mu = 1:d
        E{mu} = E{mu} * E{mu};
      end
      varargout{j}{1} = tucker (V, E);
    end
    T = T+(shat-1)*(p+1);
  end
  % optional output parameters
  if (nargout >= shat+1)
    varargout{shat+1} = s;
  end
  if (nargout >= shat+2)
    varargout{shat+2} = q;
  end
  if (nargout >= shat+3)
    varargout{shat+3} = T;
  end
%!error % in 1d, the matrix has to be A{1}, A is not allowed
%! A = randn(4);
%! phiks(0.1,A,rand(4,1),2,1e-2)
%!test % p = 0, scales
%! tau = rand;
%! A{1} = randn(4);
%! V = randn(4,1);
%! p = 0;
%! [phi,phi2] = phiks(tau,A,V,p,1e-8,2); % diffphi
%! assert(phi{1},expm(tau/2*A{1})*(expm(tau/2*A{1})*V));
%! assert(phi2{1},expm(tau/2*A{1})*V);
%! U{1} = V;
%! [phi,phi2] = phiks(tau,A,U,p,1e-8,2); % lincomb
%! assert(phi,expm(tau/2*A{1})*(expm(tau/2*A{1})*V));
%! assert(phi2,expm(tau/2*A{1})*V);
%!test % p = 0
%! tau = rand;
%! A{1} = randn(4);
%! V = randn(4,1);
%! p = 0;
%! phi = phiks(tau,A,V,p,1e-8); % diffphi
%! assert(phi{1},expm(tau*A{1})*V);
%! U{1} = V;
%! phi = phiks(tau,A,U,p,1e-8); % lincomb
%! assert(phi,expm(tau*A{1})*V);
%!test % 1d first
%! tau = rand;
%! A{1} = 10*randn(4);
%! V = randn(4,1);
%! p = 2;
%! tic
%! phi = phiks(tau,A,V,p,2^-14);
%! toc
%! I = eye(4);
%! M = A{1};
%! II = eye(4);
%! res{1} = expm(tau*M)*V;
%! res{2} = (expm(tau*M)-II)/(tau*M)*V;
%! res{3} = (expm(tau*M)-II-tau*M)/(tau*M)/(tau*M)*V;
%! assert(res,phi,-1e-3)
%! U{1} = randn(4,1);
%! U{2} = randn(4,1);
%! U{3} = randn(4,1);
%! tic
%! phicomb = phiks(tau,A,U,p,2^-14);
%! toc
%! res = (expm(tau*M)-II-tau*M)/(tau*M)/(tau*M)*U{3}+...
%!       (expm(tau*M)-II)/(tau*M)*U{2}+...
%!       expm(tau*M)*U{1};
%! assert(res,phicomb,-1e-3)
%!test % 2d
%! tau = rand;
%! A{1} = 10*randn(4);
%! A{2} = 10*randn(4);
%! V = randn(4);
%! p = 2;
%! tic
%! phi = phiks(tau,A,V,p,2^-14);
%! toc
%! I = eye(4);
%! M = kron(A{2},I)+kron(I,A{1});
%! II = eye(16);
%! res{1} = reshape(expm(tau*M)*V(:),4,4);
%! res{2} = reshape((expm(tau*M)-II)/(tau*M)*V(:),4,4);
%! res{3} = reshape((expm(tau*M)-II-tau*M)/(tau*M)/(tau*M)*V(:),4,4);
%! assert(res,phi,-1e-3)
%! U{1} = randn(4,4);
%! U{2} = randn(4,4);
%! U{3} = randn(4,4);
%! tic
%! phicomb = phiks(tau,A,U,p,2^-14);
%! toc
%! res = reshape((expm(tau*M)-II-tau*M)/(tau*M)/(tau*M)*U{3}(:)+...
%!       (expm(tau*M)-II)/(tau*M)*U{2}(:)+...
%!       expm(tau*M)*U{1}(:),4,4);
%! assert(res,phicomb,-1e-3)
%!test % 2d no shift
%! tau = rand;
%! A{1} = 10*randn(4);
%! A{2} = 10*randn(4);
%! V = randn(4);
%! p = 2;
%! tic
%! phi = phiks(tau,A,V,p,2^-14,1,0);
%! toc
%! I = eye(4);
%! M = kron(A{2},I)+kron(I,A{1});
%! II = eye(16);
%! res{1} = reshape(expm(tau*M)*V(:),4,4);
%! res{2} = reshape((expm(tau*M)-II)/(tau*M)*V(:),4,4);
%! res{3} = reshape((expm(tau*M)-II-tau*M)/(tau*M)/(tau*M)*V(:),4,4);
%! assert(res,phi,-1e-3)
%! U{1} = randn(4,4);
%! U{2} = randn(4,4);
%! U{3} = randn(4,4);
%! tic
%! phicomb = phiks(tau,A,U,p,2^-14,1,0);
%! toc
%! res = reshape((expm(tau*M)-II-tau*M)/(tau*M)/(tau*M)*U{3}(:)+...
%!       (expm(tau*M)-II)/(tau*M)*U{2}(:)+...
%!       expm(tau*M)*U{1}(:),4,4);
%! assert(res,phicomb,-1e-3)
%!test % 1d
%! tau = rand;
%! A{1} = randn(4);
%! p = 2;
%! I = eye(4);
%! M = A{1};
%! II = eye(4);
%! U{1} = randn(4,1);
%! U{2} = randn(4,1);
%! U{3} = randn(4,1);
%! tic
%! [phicomb1,phicomb2] = phiks(tau,A,U,p,2^-14,2);
%! toc
%! res1 = (expm(tau*M)-II-tau*M)/(tau*M)/(tau*M)*U{3}+...
%!       (expm(tau*M)-II)/(tau*M)*U{2}+...
%!       expm(tau*M)*U{1};
%! assert(res1,phicomb1,-1e-3)
%! tau = tau/2;
%! res2 = (expm(tau*M)-II-tau*M)/(tau*M)/(tau*M)*U{3}/4+...
%!       (expm(tau*M)-II)/(tau*M)*U{2}/2+...
%!       expm(tau*M)*U{1};
%! assert(res2,phicomb2,-1e-3)
%!test % 1d % although small norm, two output arguments required
%! tau = 0.1;
%! A{1} = randn(4);
%! p = 2;
%! I = eye(4);
%! M = A{1};
%! II = eye(4);
%! U{1} = rand(4,1)/10;
%! U{2} = rand(4,1)/10;
%! U{3} = rand(4,1)/10;
%! tic
%! [phicomb1,phicomb2] = phiks(tau,A,U,p,2^-14,2);
%! toc
%! res1 = (expm(tau*M)-II-tau*M)/(tau*M)/(tau*M)*U{3}+...
%!       (expm(tau*M)-II)/(tau*M)*U{2}+...
%!       expm(tau*M)*U{1};
%! assert(res1,phicomb1,-1e-3)
%! tau = tau/2;
%! res2 = (expm(tau*M)-II-tau*M)/(tau*M)/(tau*M)*U{3}/4+...
%!       (expm(tau*M)-II)/(tau*M)*U{2}/2+...
%!       expm(tau*M)*U{1};
%! assert(res2,phicomb2,-1e-3)
%!test % 1d
%! tau = 0.1;
%! A{1} = 10*toeplitz([2,-1,0,0]);
%! V = (1:4)';
%! p = 3;
%! tic
%! [phi,phihalf,phiquarter] = phiks(tau,A,V,p,2^-20,3);
%! toc
%! I = eye(4);
%! M = A{1};
%! II = eye(4);
%! res{1} = expm(tau*M)*V;
%! res{2} = (expm(tau*M)-II)/(tau*M)*V;
%! res{3} = (expm(tau*M)-II-tau*M)/(tau*M)/(tau*M)*V;
%! res{4} = (expm(tau*M)-II-tau*M-(tau*M)^2/2)/(tau*M)/(tau*M)/(tau*M)*V;
%! assert(res,phi,-1e-3)
%! tau = tau/2;
%! res{1} = expm(tau*M)*V;
%! res{2} = (expm(tau*M)-II)/(tau*M)*V;
%! res{3} = (expm(tau*M)-II-tau*M)/(tau*M)/(tau*M)*V;
%! res{4} = (expm(tau*M)-II-tau*M-(tau*M)^2/2)/(tau*M)/(tau*M)/(tau*M)*V;
%! assert(res,phihalf,-1e-3)
%! tau = tau/2;
%! res{1} = expm(tau*M)*V;
%! res{2} = (expm(tau*M)-II)/(tau*M)*V;
%! res{3} = (expm(tau*M)-II-tau*M)/(tau*M)/(tau*M)*V;
%! res{4} = (expm(tau*M)-II-tau*M-(tau*M)^2/2)/(tau*M)/(tau*M)/(tau*M)*V;
%! assert(res,phiquarter,-1e-3)
%!test % 1d % although small norm, two output arguments required
%! tau = 0.1;
%! A{1} = randn(4);
%! p = 2;
%! I = eye(4);
%! M = A{1};
%! II = eye(4);
%! U{1} = rand(4,1)/10;
%! U{2} = rand(4,1)/10;
%! U{3} = rand(4,1)/10;
%! tic
%! [phicomb1,phicomb2,phicomb3] = phiks(tau,A,U,p,2^-14,3);
%! toc
%! res1 = (expm(tau*M)-II-tau*M)/(tau*M)/(tau*M)*U{3}+...
%!       (expm(tau*M)-II)/(tau*M)*U{2}+...
%!       expm(tau*M)*U{1};
%! assert(res1,phicomb1,-1e-3)
%! tau = tau/2;
%! res2 = (expm(tau*M)-II-tau*M)/(tau*M)/(tau*M)*U{3}/4+...
%!       (expm(tau*M)-II)/(tau*M)*U{2}/2+...
%!       expm(tau*M)*U{1};
%! assert(res2,phicomb2,-1e-3)
%! tau = tau/2;
%! res3 = (expm(tau*M)-II-tau*M)/(tau*M)/(tau*M)*U{3}/16+...
%!       (expm(tau*M)-II)/(tau*M)*U{2}/4+...
%!       expm(tau*M)*U{1};
%! assert(res3,phicomb3,-1e-3)
%!test % 1d first
%! tau = rand;
%! A{1} = 10*randn(4);
%! V = randn(4,1);
%! p = 3;
%! tic
%! phi = phiks(tau,A,V,p,2^-14);
%! toc
%! I = eye(4);
%! M = A{1};
%! II = eye(4);
%! res{1} = expm(tau*M)*V;
%! res{2} = (expm(tau*M)-II)/(tau*M)*V;
%! res{3} = (expm(tau*M)-II-tau*M)/(tau*M)/(tau*M)*V;
%! res{4} = (expm(tau*M)-II-tau*M-(tau*M)^2/2)/(tau*M)/(tau*M)/(tau*M)*V;
%! assert(res,phi,-1e-3)
%! U{1} = randn(4,1);
%! U{2} = randn(4,1);
%! U{3} = randn(4,1);
%! U{4} = randn(4,1);
%! tic
%! phicomb = phiks(tau,A,U,p,2^-14);
%! toc
%! res = (expm(tau*M)-II-tau*M-(tau*M)^2/2)/(tau*M)/(tau*M)/(tau*M)*U{4}+...
%!       (expm(tau*M)-II-tau*M)/(tau*M)/(tau*M)*U{3}+...
%!       (expm(tau*M)-II)/(tau*M)*U{2}+...
%!       expm(tau*M)*U{1};
%! assert(res,phicomb,-1e-3)
%!test % 2d complex
%! tau = rand;
%! A{1} = randn(4)+1i*10*randn(4);
%! A{2} = randn(4)+1i*10*randn(4);
%! V = randn(4)+1i*randn(4);
%! p = 2;
%! tic
%! phi = phiks(tau,A,V,p,2^-14);
%! toc
%! I = eye(4);
%! M = kron(A{2},I)+kron(I,A{1});
%! II = eye(16);
%! res{1} = reshape(expm(tau*M)*V(:),4,4);
%! res{2} = reshape((expm(tau*M)-II)/(tau*M)*V(:),4,4);
%! res{3} = reshape((expm(tau*M)-II-tau*M)/(tau*M)/(tau*M)*V(:),4,4);
%! assert(res,phi,-1e-3)
%! U{1} = randn(4)+1i*randn(4,4);
%! U{2} = randn(4)+1i*randn(4,4);
%! U{3} = randn(4)+1i*randn(4,4);
%! tic
%! phicomb = phiks(tau,A,U,p,2^-14);
%! toc
%! res = reshape((expm(tau*M)-II-tau*M)/(tau*M)/(tau*M)*U{3}(:)+...
%!       (expm(tau*M)-II)/(tau*M)*U{2}(:)+...
%!       expm(tau*M)*U{1}(:),4,4);
%! assert(res,phicomb,-1e-3)
%!test % 1d % although small norm, four output arguments required
%! tau = 0.1;
%! A{1} = magic(4)+eye(4);
%! p = 2;
%! I = eye(4);
%! M = A{1};
%! II = eye(4);
%! U{1} = [1:4]'/70;
%! U{2} = [2:5]'/70;
%! U{3} = [3:6]'/70;
%! tic
%! [phicomb1,phicomb2,s,q] = phiks(tau,A,U,p,2^-14,2);
%! toc
%! res1 = (expm(tau*M)-II-tau*M)/(tau*M)/(tau*M)*U{3}+...
%!       (expm(tau*M)-II)/(tau*M)*U{2}+...
%!       expm(tau*M)*U{1};
%! assert(res1,phicomb1,-1e-3)
%! tau = tau/2;
%! res2 = (expm(tau*M)-II-tau*M)/(tau*M)/(tau*M)*U{3}/4+...
%!       (expm(tau*M)-II)/(tau*M)*U{2}/2+...
%!       expm(tau*M)*U{1};
%! assert(res2,phicomb2,-1e-3)
%! assert(s,2)
%! assert(q,3)
%!test % 1d
%! tau = 0.1;
%! A{1} = 10*toeplitz([2,-1,0,0]);
%! V = (1:4)';
%! p = 3;
%! tic
%! [phi,phihalf,phiquarter,s,q] = phiks(tau,A,V,p,2^-20,3);
%! toc
%! I = eye(4);
%! M = A{1};
%! II = eye(4);
%! res{1} = expm(tau*M)*V;
%! res{2} = (expm(tau*M)-II)/(tau*M)*V;
%! res{3} = (expm(tau*M)-II-tau*M)/(tau*M)/(tau*M)*V;
%! res{4} = (expm(tau*M)-II-tau*M-(tau*M)^2/2)/(tau*M)/(tau*M)/(tau*M)*V;
%! assert(res,phi,-1e-3)
%! tau = tau/2;
%! res{1} = expm(tau*M)*V;
%! res{2} = (expm(tau*M)-II)/(tau*M)*V;
%! res{3} = (expm(tau*M)-II-tau*M)/(tau*M)/(tau*M)*V;
%! res{4} = (expm(tau*M)-II-tau*M-(tau*M)^2/2)/(tau*M)/(tau*M)/(tau*M)*V;
%! assert(res,phihalf,-1e-3)
%! tau = tau/2;
%! res{1} = expm(tau*M)*V;
%! res{2} = (expm(tau*M)-II)/(tau*M)*V;
%! res{3} = (expm(tau*M)-II-tau*M)/(tau*M)/(tau*M)*V;
%! res{4} = (expm(tau*M)-II-tau*M-(tau*M)^2/2)/(tau*M)/(tau*M)/(tau*M)*V;
%! assert(res,phiquarter,-1e-3)
%! assert(s,2);
%! assert(q,5);
%!test % 1d
%! tau = 0.1;
%! A{1} = 200*toeplitz([2,-1,0,0]);
%! V = (1:4)';
%! p = 3;
%! tic
%! [phi,phihalf,phiquarter,s,q] = phiks(tau,A,V,p,2^-20,3);
%! toc
%! I = eye(4);
%! M = A{1};
%! II = eye(4);
%! res{1} = expm(tau*M)*V;
%! res{2} = (expm(tau*M)-II)/(tau*M)*V;
%! res{3} = (expm(tau*M)-II-tau*M)/(tau*M)/(tau*M)*V;
%! res{4} = (expm(tau*M)-II-tau*M-(tau*M)^2/2)/(tau*M)/(tau*M)/(tau*M)*V;
%! assert(res,phi,-1e-3)
%! tau = tau/2;
%! res{1} = expm(tau*M)*V;
%! res{2} = (expm(tau*M)-II)/(tau*M)*V;
%! res{3} = (expm(tau*M)-II-tau*M)/(tau*M)/(tau*M)*V;
%! res{4} = (expm(tau*M)-II-tau*M-(tau*M)^2/2)/(tau*M)/(tau*M)/(tau*M)*V;
%! assert(res,phihalf,-1e-3)
%! tau = tau/2;
%! res{1} = expm(tau*M)*V;
%! res{2} = (expm(tau*M)-II)/(tau*M)*V;
%! res{3} = (expm(tau*M)-II-tau*M)/(tau*M)/(tau*M)*V;
%! res{4} = (expm(tau*M)-II-tau*M-(tau*M)^2/2)/(tau*M)/(tau*M)/(tau*M)*V;
%! assert(res,phiquarter,-1e-3)
%! assert(s,4);
%! assert(q,7);
%!test % 2d
%! p =2;
%! tau = 1/2;
%! A{1} = (magic(4)+eye(4))/100;
%! A{2} = (magic(4)+eye(4))/100;
%! M = kron(eye(4),A{1})+kron(A{2},eye(4));
%! II = eye(16);
%! U{1} = toeplitz([-2,1,0,0]);
%! U{2} = toeplitz([-2,1,0,0]);
%! U{3} = toeplitz([-2,1,0,0]);
%! tic
%! [phicomb,s] = phiks(tau,A,U,p,2^-14);
%! toc
%! res = reshape((expm(tau*M)-II-tau*M)/(tau*M)/(tau*M)*U{3}(:)+...
%!       (expm(tau*M)-II)/(tau*M)*U{2}(:)+...
%!       expm(tau*M)*U{1}(:),4,4);
%! assert(res,phicomb,-1e-3)
%! assert(s,0)
%!test % V{2} = 0
%! A{1} = randn(2);
%! A{2} = randn(2);
%! A{3} = randn(2);
%! V{1} = randn(2,2,2);
%! V{2} = zeros(2,2,2);
%! V{3} = randn(2,2,2);
%! [ret,s,q,c] = phiks(1,A,V,2,1e-8);
%! V{2} = 0;
%! [ret0,s0,q0,c0] = phiks(1,A,V,2,1e-8);
%! assert(ret,ret0);
%! assert(c,c0+1);
%!test % V{1} = 0
%! A{1} = randn(2);
%! A{2} = randn(2);
%! A{3} = randn(2);
%! V{1} = zeros(2,2,2);
%! V{2} = randn(2,2,2);
%! V{3} = randn(2,2,2);
%! [ret,s,q,c] = phiks(1,A,V,2,1e-8);
%! V{1} = 0;
%! [ret0,s0,q0,c0] = phiks(1,A,V,2,1e-8);
%! assert(ret,ret0);
%! assert(c,c0+1);
%!test % V{1} = V{2} = 0 and scales>1
%! A{1} = 10*randn(2);
%! A{2} = 10*randn(2);
%! A{3} = 10*randn(2);
%! V{1} = zeros(2,2,2);
%! V{2} = zeros(2,2,2);
%! V{3} = randn(2,2,2);
%! V{4} = randn(2,2,2);
%! [ret,rethalf,s,q,c] = phiks(1,A,V,3,1e-8,2);
%! V{1} = 0;
%! V{2} = 0;
%! [ret0,rethalf0,s0,q0,c0] = phiks(1,A,V,3,1e-8,2);
%! assert(ret,ret0);
%! assert(rethalf,rethalf0);
%! assert(c,c0+3);
%!test % number of Tucker
%! A{1} = toeplitz([-2,1,0,0]);
%! A{2} = toeplitz([-2,1,0,0]);
%! V = magic(4);
%! [res,s,q,T] = phiks(1,A,V,0,1e-8);
%! assert(T,1)
%! [res1,res2,s,q,T] = phiks(1,A,V,0,1e-8,2);
%! assert(T,2)
%! [res,s,q,T] = phiks(1,A,V,1,1e-8);
%! assert(T,8)
%! [res1,res2,s,q,T] = phiks(1,A,V,1,1e-8,2);
%! assert(T,9)
%! [res1,res2,s,q,T] = phiks(10,A,V,1,1e-8,2);
%! assert(T,12)
%! clear V;
%! V{1} = magic(4);
%! V{2} = magic(4);
%! V{3} = magic(4);
%! [res,s,q,T] = phiks(1,A,V,0,1e-8);
%! assert(T,1);
%! [res1,res2,s,q,T] = phiks(1,A,V,0,1e-8,2);
%! assert(T,2);
%! [res,s,q,T] = phiks(1,A,V,1,1e-8);
%! assert(T,8);
%! [res1,res2,s,q,T] = phiks(1,A,V,1,1e-8,2);
%! assert(T,9);
%! [res,s,q,T] = phiks(1,A,V,2,1e-8);
%! assert(T,14);
%! [res1,res2,s,q,T] = phiks(1,A,V,2,1e-8,2);
%! assert(T,15);
%! V{1} = 0;
%! [res1,res2,s,q,T] = phiks(1,A,V,2,1e-8,2);
%! assert(T,13);
%! V{2} = 0;
%! [res1,res2,s,q,T] = phiks(1,A,V,2,1e-8,2);
%! assert(T,12);
%!test % 3d, complex tau
%! tau = 0.1+1i*0.1;
%! A{1} = 20*toeplitz([2,-1,0,0]);
%! A{2} = 20*toeplitz([2,-1,0,0]);
%! A{3} = 20*toeplitz([2,-1,0,0]);
%! V = reshape(1:64,4,4,4);
%! p = 3;
%! tic
%! [phi,phihalf,phiquarter,s,q] = phiks(tau,A,V,p,2^-20,3);
%! toc
%! I = eye(4);
%! M = kron(I,kron(I,A{1}))+kron(I,kron(A{2},I))+kron(A{3},kron(I,I));
%! II = eye(64);
%! res{1} = expm(tau*M)*V(:);
%! res{2} = (expm(tau*M)-II)/(tau*M)*V(:);
%! res{3} = (expm(tau*M)-II-tau*M)/(tau*M)/(tau*M)*V(:);
%! res{4} = (expm(tau*M)-II-tau*M-(tau*M)^2/2)/(tau*M)/(tau*M)/(tau*M)*V(:);
%! assert(res,{phi{1}(:),phi{2}(:),phi{3}(:),phi{4}(:)},-1e-3)
%! tau = tau/2;
%! res{1} = expm(tau*M)*V(:);
%! res{2} = (expm(tau*M)-II)/(tau*M)*V(:);
%! res{3} = (expm(tau*M)-II-tau*M)/(tau*M)/(tau*M)*V(:);
%! res{4} = (expm(tau*M)-II-tau*M-(tau*M)^2/2)/(tau*M)/(tau*M)/(tau*M)*V(:);
%! assert(res,{phihalf{1}(:),phihalf{2}(:),phihalf{3}(:),phihalf{4}(:)},-1e-3)
%! tau = tau/2;
%! res{1} = expm(tau*M)*V(:);
%! res{2} = (expm(tau*M)-II)/(tau*M)*V(:);
%! res{3} = (expm(tau*M)-II-tau*M)/(tau*M)/(tau*M)*V(:);
%! res{4} = (expm(tau*M)-II-tau*M-(tau*M)^2/2)/(tau*M)/(tau*M)/(tau*M)*V(:);
%! assert(res,{phiquarter{1}(:),phiquarter{2}(:),...
%!        phiquarter{3}(:),phiquarter{4}(:)},-1e-3)
%! assert(s,2);
%! assert(q,10);
%!demo % large 3d
%! n = 150;
%! A{1} = toeplitz([-2,1,zeros(1,n-2)]);
%! A{2} = A{1};
%! A{3} = A{1};
%! V{1} = ones(n,n,n);
%! V{2} = V{1};
%! V{3} = V{1};
%! V{4} = V{1};
%! tic
%! [ret,rethalf,s,q,c] = phiks(1,A,V,3,1e-8,2);
%! toc
%! V = ones(n,n,n);
%! tic
%! [ret,rethalf,s,q,c] = phiks(1,A,V,3,1e-8,2);
%! toc
