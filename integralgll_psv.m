function [I, E, T] = integralgll_psv (tau, A, V, p, s, q, sigma)
%
% function [I, E] = integralgll_psv (tau, A, V, p, s, q, sigma)
%
% Gauss-Lobatto-Legendre quadrature formula with q nodes, for the
% computation of
%
%  I{ell+1} = phi_ell(tau*K/2^s)*V, ell = 0,1,2,...,p
%
% where K is the Kronecker sum of A{d}+sigma(d)*I,...,A{1}+sigma(d)*I.
% Moreover, E{mu} = exp(tau*A{mu}/2^s)*exp(tau*sigma(mu)/2^s).
[theta,w] = XW(q);
tau2s = tau / 2 ^ s;
d = length (A);
i = q; % theta(i) = 1
weight = w(i);
I{2} = V*weight;
for ell = 2:p
  weight = weight*theta(i)/(ell-1);
  I{ell+1} = V*weight;
end
i = q-1;
for mu = 1:d
  E{mu} = expm ((1 - theta(i)) * tau2s * A{mu}) * ...
          exp ((1 - theta(i)) * tau2s * sigma(mu));
end
B = tucker (V, E);
weight = w(i);
I{2} = I{2}+B*weight;
for ell = 2:p
  weight = weight*theta(i)/(ell-1);
  I{ell+1} = I{ell+1}+B*weight;
end
for i = q-2:-1:2
  for mu = 1:d
    E{mu} = E{mu} * expm ((theta(i+1)-theta(i)) * tau2s * A{mu}) * ...
            exp ((theta(i+1) - theta(i)) * tau2s * sigma(mu));
  end
  B = tucker (V, E);
  weight = w(i);
  I{2} = I{2}+B*weight;
  for ell = 2:p
    weight = weight*theta(i)/(ell-1);
    I{ell+1} = I{ell+1}+B*weight;
  end
end
i = 1; % theta(i) = 0
for mu = 1:d
  E{mu} = E{mu} * expm ((theta(i+1)-theta(i)) * tau2s * A{mu}) * ...
          exp ((theta(i+1) - theta(i)) * tau2s * sigma(mu));
end
I{1} = tucker (V, E);
I{2} = I{2}+I{1}*w(i);
T = q-1;
%!test
%! A{1} = 0*randn(4)+0*1i*randn(4);
%! A{2} = 0*randn(4)+0*1i*randn(4);
%! K = kron(eye(4),A{1})+kron(A{2},eye(4));
%! tau = rand;
%! p = 4;
%! V = randn(4,4);
%! s = 3;
%! q = 3;
%! sigma = [0,0];
%! [I,E] = integralgll_psv(tau, A, V, p, s, q, sigma);
%! assert(I{1}(:),expm(tau*K/2^s)*V(:),1e-14)
%! assert(I{2}(:),phi1m(tau*K/2^s)*V(:),1e-14)
%! assert(I{3}(:),phi2m(tau*K/2^s)*V(:),1e-14)
%! assert(I{4}(:),phi3m(tau*K/2^s)*V(:),1e-14)
%! assert(I{5}(:),phi4m(tau*K/2^s)*V(:),1e-14)
%! assert(E{1},expm(tau*A{1}/2^s),1e-14)
%! assert(E{2},expm(tau*A{2}/2^s),1e-14)
%!test
%! A{1} = randn(4)+1i*randn(4);
%! A{2} = randn(4)+1i*randn(4);
%! K = kron(eye(4),A{1})+kron(A{2},eye(4));
%! tau = rand;
%! p = 4;
%! V = randn(4,4);
%! s = 3;
%! q = 12;
%! sigma = [0,0];
%! [I,E] = integralgll_psv(tau, A, V, p, s, q, sigma);
%! assert(I{1}(:),expm(tau*K/2^s)*V(:),1e-14)
%! assert(I{2}(:),phi1m(tau*K/2^s)*V(:),1e-14)
%! assert(I{3}(:),phi2m(tau*K/2^s)*V(:),1e-14)
%! assert(I{4}(:),phi3m(tau*K/2^s)*V(:),1e-14)
%! assert(I{5}(:),phi4m(tau*K/2^s)*V(:),1e-14)
%! assert(E{1},expm(tau*A{1}/2^s),1e-14)
%! assert(E{2},expm(tau*A{2}/2^s),1e-14)
%!test
%! A{1} = randn(4)+1i*randn(4);
%! A{2} = randn(4)+1i*randn(4);
%! K = kron(eye(4),A{1})+kron(A{2},eye(4));
%! sigma(1) = trace(A{1})/4;
%! A{1} = A{1}-sigma(1)*eye(4);
%! sigma(2) = trace(A{2})/4;
%! A{2} = A{2}-sigma(2)*eye(4);
%! tau = rand;
%! p = 4;
%! V = randn(4,4);
%! s = 3;
%! q = 12;
%! [I,E] = integralgll_psv(tau, A, V, p, s, q, sigma);
%! assert(I{1}(:),expm(tau*K/2^s)*V(:),1e-14)
%! assert(I{2}(:),phi1m(tau*K/2^s)*V(:),1e-14)
%! assert(I{3}(:),phi2m(tau*K/2^s)*V(:),1e-14)
%! assert(I{4}(:),phi3m(tau*K/2^s)*V(:),1e-14)
%! assert(I{5}(:),phi4m(tau*K/2^s)*V(:),1e-14)
%! assert(E{1},expm(tau*A{1}/2^s)*exp(tau*sigma(1)/2^s),1e-14)
%! assert(E{2},expm(tau*A{2}/2^s)*exp(tau*sigma(2)/2^s),1e-14)
