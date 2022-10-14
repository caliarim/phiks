function K = kronsum(varargin)
% KRONSUM Kronecker sum of matrices.
%    K = KRONSUM(L) produces the sparse matrix
%
%    KRON(I, KRON(I, ..., KRON(I, L{1})))+KRON(I, KRON(I, ..., KRON(L{2}, I)))+
%        ... KRON(L{d}, KRON(I, ..., KRON(I, I)))
%
%    where the cell array L contains complex matrices L{mu} of size
%    n_{mu} x n_{mu}, while I are identity matrices of suitable sizes. In other
%    words, we have
%
%    K = L{d} \oplus ... \oplus L{1},
%
%    with d > 1.
%
%    K = KRONSUM(L1, L2, ..., Ld) produces the Kronecker sum of the
%    complex matrices L1, L2, ..., Ld, with Lmu of size n_{mu} x n_{mu}.
%
%    [CCZ22a] M. Caliari, F. Cassini, and F. Zivcovich,
%            A mu-mode BLAS approach for multidimensional tensor-structured
%            problems, Numerical Algorithms, 2022
  if (nargin < 1)
    error('Not enough input arguments.')
  end
  if (iscell(varargin{1}))
    varargin = varargin{1};
  end
  d = length(varargin);
  if (any(cellfun(@isempty,varargin)) || (d == 1))
    error('Invalid input arguments.')
  end
  murange = 1:length(varargin);
  murange = murange(~cellfun(@isempty,varargin));
  n = cellfun(@length,varargin);
  K = kron(speye(prod(n(2:d))), varargin{1});
  for mu = 2:d-1
    K = K+kron(kron(speye(prod(n(mu+1:d))), varargin{mu}), ...
               speye(prod(n(1:mu-1))));
  end
  K = K+kron(varargin{d}, speye(prod(n(1:d-1))));
end
%!test % different input form
%! A{1} = sprand(3,3,0.1);
%! A{2} = sprand(4,4,0.1);
%! assert(kronsum(A),kronsum(A{1},A{2}))
%!test % 2d
%! A{1} = sprand(3,3,0.1);
%! A{2} = sprand(4,4,0.1);
%! assert(kronsum(A),kron(speye(4),A{1})+kron(A{2},speye(3)))
%!test % 2d full
%! A{1} = rand(3,3);
%! A{2} = rand(4,4);
%! assert(kronsum(A),kron(speye(4),A{1})+kron(A{2},speye(3)))
%!test % 3d
%! A{1} = sprand(3,3,0.1);
%! A{2} = sprand(4,4,0.1);
%! A{3} = sprand(5,5,0.1);
%! assert(kronsum(A),...
%!  kron(kron(speye(5),speye(4)),A{1})+...
%!  kron(kron(speye(5),A{2}),speye(3))+...
%!  kron(kron(A{3},speye(4)),speye(3)),eps)
%!test % 3d full
%! A{1} = rand(3,3);
%! A{2} = rand(4,4);
%! A{3} = rand(5,5);
%! assert(kronsum(A),...
%!  kron(kron(speye(5),speye(4)),A{1})+...
%!  kron(kron(speye(5),A{2}),speye(3))+...
%!  kron(kron(A{3},speye(4)),speye(3)),2*eps)
%! assert(issparse(kronsum(A)),true)
%!test %kronsumv
%! A = sprand(3,3,0.1);
%! B = sprand(4,4,0.1);
%! C = sprand(5,5,0.1);
%! U = randn(3,4,5);
%! assert(kronsum(A,B,C)*U(:), kronsumv(U, A, B, C)(:),10*eps)
%!test % complex
%! A = randn(3)+1i*randn(3);
%! B = randn(2)+1i*randn(2);
%! ref = kron(speye(2),A)+kron(B,speye(3));
%! assert(kronsum(A,B), ref)
%!error
%! kronsum()
%!error
%! kronsum({rand(2),[]})
%!error
%! kronsum(rand(3));
