function T = mump(T, L, mu)
% MUMP mu-mode product.
%     S = MUMP(T, L, mu) computes the mu-mode product of the complex tensor T
%     of size m_1 x ... x m_d with the complex matrix L of size n_{mu} x m_{mu},
%     that is
%
%     S = T x_{mu} L
%
%    [CCZ22a] M. Caliari, F. Cassini, and F. Zivcovich,
%            A mu-mode BLAS approach for multidimensional tensor-structured
%            problems, Numerical Algorithms, 2022
  if (nargin < 3)
    error('Not enough input arguments.');
  end
  if (isempty(T) || isempty(L) || isempty(mu))
    error('Not enough non-empty input arguments');
  end
  sT = [size(T), ones(1, mu-length(size(T)))];
  sL = size(L);
  sT(mu) = sL(1);
  if (mu == 1)
    T = reshape(L*reshape(T, sL(2), []), sT);
  elseif (mu == length(sT))
    T = reshape(reshape(T, [], sL(2))*L.', sT);
  else
    idx = [mu, 1:(mu-1), (mu+1):length(sT)];
    T = ipermute(reshape(L*reshape(permute(T, idx), sL(2), []), sT(idx)), idx);
  end
end
%!test % 1d
%! T = (1:4)';
%! A = reshape (1:8,2,4);
%! Y = mump(T,A,1);
%! assert(Y,A*T)
%!test %2d
%! T = reshape(1:16,4,4);
%! A = reshape (1:8,2,4);
%! Y1 = mump(T,A,1);
%! assert(Y1,A*T)
%! Y2 = mump(T,A,2);
%! assert(Y2,(A*T')')
%!test %2d rect
%! T = randn(2,3);
%! A = randn(4,2);
%! B = randn(5,3);
%! Y1 = mump(T,A,1);
%! assert(Y1,A*T)
%! Y2 = mump(T,B,2);
%! assert(Y2,(B*T')')
%!test %2d rect complex
%! T = randn(2,3)+1i*randn(2,3);
%! A = randn(4,2)+1i*randn(4,2);
%! B = randn(5,3)+1i*randn(5,3);
%! Y1 = mump(T,A,1);
%! assert(Y1,A*T)
%! Y2 = mump(T,B,2);
%! assert(Y2,(B*T.').')
%!test % 3d
%! T(:,:,1) = [1,2,3;4,5,6];
%! T(:,:,2) = [7,8,9;10,11,12];
%! T(:,:,3) = [13,14,15;16,17,18];
%! T(:,:,4) = [19,20,21;22,23,24];
%! A1 = reshape(1:4,2,2);
%! Y1 = mump(T,A1,1);
%! ref(:,:,1) = [13,17,21;18,24,30];
%! ref(:,:,2) = [37,41,45;54,60,66];
%! ref(:,:,3) = [61,65,69;90,96,102];
%! ref(:,:,4) = [85,89,93;126,132,138];
%! assert(Y1,ref)
%! A2 = reshape(1:9,3,3);
%! Y2 = mump(T,A2,2);
%! ref(:,:,1) = [30,36,42;66,81,96];
%! ref(:,:,2) = [102,126,150;138,171,204];
%! ref(:,:,3) = [174,216,258;210,261,312];
%! ref(:,:,4) = [246,306,366;282,351,420];
%! assert(Y2,ref)
%!test
%! A{1} = randn(2);
%! I{1} = eye(2);
%! A{2} = randn(3);
%! I{2} = eye(3);
%! A{3} = randn(4);
%! I{3} = eye(4);
%! AA = kron(I{3},kron(I{2},A{1}))+...
%!      kron(I{3},kron(A{2},I{1}))+...
%!      kron(A{3},kron(I{2},I{1}));
%! U0 = randn(2,3,4);
%! Y = mump(U0,A{1},1)+...
%!     mump(U0,A{2},2)+...
%!     mump(U0,A{3},3);
%! assert(AA * U0(:),Y(:),1e-10)
%!error
%! mump();
%!error
%! mump(rand(2,3,4))
%!error
%! mump(rand(2,3,4),rand(2))
%!error
%! mump([],rand(2),1)
%!error
%! mump(randn(2),[],1)
%!error
%! mump(randn(2),rand(2),[])
