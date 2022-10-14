function [I, E, T] = integralgll_lcp (tau, A, V, p, s, q, sigma, zeroV)
%
% function [I,E] = integralgll_lcp (tau, A, V, p, s, q, sigma)
%
% Gauss-Lobatto-Legendre quadrature formula with q nodes, for the
% computation of
%
%  I{ell} = Phi_s(tau*K)V{p-ell+2:p}, ell = 1,2,...,p
%
% where K is the Kronecker sum of A{d}+sigma(d)*I,...,A{1}+sigma(d)*I.
% Moreover, E{mu} = exp(tau*A{mu}/2^s)*exp(tau*sigma(mu)/2^s).
[theta,w] = XW(q);
tau2s = tau / 2 ^ s;
w2s = w / 2 ^ s;
d = length (A);
i = q; % theta(i) = 1
for ell = 1+(s==0)*(p-1):p % if s == 0, compute only B{p}
  B = V{p+2-ell};
  fact = theta(i) / 2 ^ s;
  for k = ell-1:-1:1
    B = B+V{p+2-k}*fact;
    fact = fact * theta(i) / 2 ^ s / (ell - k + 1);
  end
  I{ell} = B*w2s(i);
end
i = q-1;
for mu = 1:d
  E{mu} = expm ((1 - theta(i)) * tau2s * A{mu}) * ...
          exp ((1 - theta(i)) * tau2s * sigma(mu));
end
for ell = 1+(s==0)*(p-1):p % if s == 0, compute only B{p}
  B = V{p+2-ell};
  fact = theta(i) / 2 ^ s;
  for k = ell-1:-1:1
    B = B+V{p+2-k}*fact;
    fact = fact * theta(i) / 2 ^ s / (ell - k + 1);
  end
  B = tucker (B, E);
  I{ell} = I{ell}+B*w2s(i);
end
for i = q-2:-1:2
  for mu = 1:d
    E{mu} = E{mu}*expm ((theta(i+1) - theta(i)) * tau2s * A{mu}) * ...
            exp ((theta(i+1) - theta(i)) * tau2s * sigma(mu));
  end
  for ell = 1+(s==0)*(p-1):p % if s == 0, compute only B{p}
    B = V{p+2-ell};
    fact = theta(i) / 2 ^ s;
    for k = ell-1:-1:1
      B = B+V{p+2-k}*fact;
      fact = fact * theta(i) / 2 ^ s / (ell - k + 1);
    end
    B = tucker (B, E);
    I{ell} = I{ell}+B*w2s(i);
  end
end
i = 1; % theta(i) = 0
for mu = 1:d
  E{mu} = E{mu}*expm ((theta(i+1) - theta(i)) * tau2s * A{mu}) * ...
          exp ((theta(i+1) - theta(i)) * tau2s * sigma(mu));
end
T = (q-2)*(1+(s~=0)*(p-1));
for ell = 1+(s==0)*(p-1):p % if s == 0, compute only B{p}
  if (~zeroV(p+2-ell))  % skip Tucker with zeros
    I{ell} = I{ell}+tucker(V{p+2-ell}, E)*w2s(i);
    T = T+1;
  end
end
%!test
%! A{1} = (randn(4)+1i*randn(4));
%! A{2} = (randn(4)+1i*randn(4));
%! K = kron(eye(4),A{1})+kron(A{2},eye(4));
%! tau = rand;
%! p = 4;
%! V{1} = randn(4,4);
%! V{2} = randn(4,4);
%! V{3} = randn(4,4);
%! V{4} = randn(4,4);
%! V{5} = randn(4,4);
%! zeroV = false(1,5);
%! s = 2;
%! q = 12;
%! sigma = [0,0];
%! [I,E] = integralgll_lcp(tau, A, V, p, s, q, sigma, zeroV);
%! assert(I{4}(:),phi1m(tau*K/2^s)*V{2}(:)/2^s+...
%!                phi2m(tau*K/2^s)*V{3}(:)/2^(2*s)+...
%!                phi3m(tau*K/2^s)*V{4}(:)/2^(3*s)+...
%!                phi4m(tau*K/2^s)*V{5}(:)/2^(4*s),1e-14)
%! assert(I{3}(:),phi1m(tau*K/2^s)*V{3}(:)/2^s+...
%!                phi2m(tau*K/2^s)*V{4}(:)/2^(2*s)+...
%!                phi3m(tau*K/2^s)*V{5}(:)/2^(3*s),1e-14)
%! assert(I{2}(:),phi1m(tau*K/2^s)*V{4}(:)/2^s+...
%!                phi2m(tau*K/2^s)*V{5}(:)/2^(2*s),1e-14)
%! assert(I{1}(:),phi1m(tau*K/2^s)*V{5}(:)/2^s,1e-14)
%! assert(E{1},expm(tau*A{1}/2^s),1e-14)
%! assert(E{2},expm(tau*A{2}/2^s),1e-14)
%!test
%! A{1} = zeros(4);
%! A{2} = zeros(4);
%! K = kron(eye(4),A{1})+kron(A{2},eye(4));
%! tau = rand;
%! p = 4;
%! V{1} = randn(4,4);
%! V{2} = randn(4,4);
%! V{3} = randn(4,4);
%! V{4} = randn(4,4);
%! V{5} = randn(4,4);
%! zeroV = false(1,5);
%! s = 1;
%! q = 3;
%! sigma = [0,0];
%! [I,E] = integralgll_lcp(tau, A, V, p, s, q, sigma, zeroV);
%! assert(I{4}(:),phi1m(tau*K/2^s)*V{2}(:)/2^s+...
%!                phi2m(tau*K/2^s)*V{3}(:)/2^(2*s)+...
%!                phi3m(tau*K/2^s)*V{4}(:)/2^(3*s)+...
%!                phi4m(tau*K/2^s)*V{5}(:)/2^(4*s),1e-14)
%! assert(I{3}(:),phi1m(tau*K/2^s)*V{3}(:)/2^s+...
%!                phi2m(tau*K/2^s)*V{4}(:)/2^(2*s)+...
%!                phi3m(tau*K/2^s)*V{5}(:)/2^(3*s),1e-14)
%! assert(I{2}(:),phi1m(tau*K/2^s)*V{4}(:)/2^s+...
%!                phi2m(tau*K/2^s)*V{5}(:)/2^(2*s),1e-14)
%! assert(I{1}(:),phi1m(tau*K/2^s)*V{5}(:)/2^s,1e-14)
%! assert(E{1},expm(tau*A{1}/2^s),1e-14)
%! assert(E{2},expm(tau*A{2}/2^s),1e-14)
