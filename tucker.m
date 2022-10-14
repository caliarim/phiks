function T = tucker(T,varargin)
% TUCKER Tucker operator.
%    S = TUCKER(T, L) computes the Tucker operator
%
%       S = T x_1 L{1} x_2 L{2} x_3 ... x_d L{d}.
%
%    Here T is a complex tensor of size m_1 x ... x m_d, L a cell array
%    of complex matrices (L{mu} of size n_{mu} x m_{mu}) and x_mu denotes
%    the mu-mode product.
%
%    S = TUCKER(T, L1, L2, ..., Ld) computes the Tucker operator
%
%       S = T x_1 L1 x_2 L2 x_3 ... x_d Ld.
%
%    Here T is a complex tensor of size m_1 x ... x m_d, while Lmu is a
%    complex matrix of size n_{mu} x m_{mu}.
%
%    In both cases, if the entry corresponding to the mu-th matrix is empty,
%    then the associated mu-mode product is skipped.
%
%    [CCZ22a] M. Caliari, F. Cassini, and F. Zivcovich,
%             A mu-mode BLAS approach for multidimensional tensor-structured
%             problems, Numerical Algorithms 2022
  if (nargin < 2)
    error('Not enough input arguments.');
  end
  if (iscell (varargin{1}))
    varargin = varargin{1};
  end
  eidx = ~cellfun(@isempty, varargin);
  sT = [size(T), ones(1, length(varargin)-find(flip(eidx), 1)+1-length(size(T)))];
  lT = length(sT);
  mur = 1:length(varargin);
  mur = mur(eidx);
  lmu = length(mur);
  if (lmu == 0)
    error('Not enough non-empty input arguments.');
  end
  if (mur(1) == 1) && (mur(lmu) == lT)
    T = varargin{1}*reshape(T, sT(1), []);
    sT(1) = size(T, 1);
    T = reshape(T, [], sT(mur(lmu)))*varargin{mur(lmu)}.';
    sT(mur(lmu)) = size(T, 2);
    T = reshape(T, sT);
    mur = mur(2:lmu-1);
    lmu = lmu-2;
  elseif (mur(1) == 1) && (mur(lmu) ~= lT)
    T = varargin{1}*reshape(T, sT(1), []);
    sT(1) = size(T, 1);
    T = reshape(T, sT);
    mur = mur(2:lmu);
    lmu = lmu-1;
  elseif (mur(1) ~= 1) && (mur(lmu) == lT)
    T = reshape(T, [], sT(mur(lmu)))*varargin{mur(lmu)}.';
    sT(mur(lmu)) = size(T, 2);
    T = reshape(T, sT);
    mur = mur(1:lmu-1);
    lmu = lmu-1;
  end
  if (lmu > 0)
    T = permute(T, [mur(1), 1:(mur(1)-1), (mur(1)+1):lT]);
    for mu = 1:(lmu-1)
      T = varargin{mur(mu)}*reshape(T, sT(mur(mu)), []);
      sT(mur(mu)) = size(T, 1);
      T = permute(reshape(T, sT([mur(mu), 1:(mur(mu)-1), (mur(mu)+1):lT])), ...
      [mur(mu+1), 2:mur(mu), 1, (mur(mu)+1):(mur(mu+1)-1), (mur(mu+1)+1):lT]);
    end
    T = varargin{mur(lmu)}*reshape(T, sT(mur(lmu)), []);
    sT(mur(lmu)) = size(T, 1);
    T = ipermute(reshape(T, sT([mur(lmu), 1:(mur(lmu)-1), (mur(lmu)+1):lT])), ...
        [mur(lmu), 1:(mur(lmu)-1), (mur(lmu)+1):lT]);
  end
end
%!test % different input form
%! T = randn(2,3,4);
%! A{1} = randn(2);
%! A{2} = randn(3);
%! A{3} = randn(4);
%! assert(tucker(T,A),tucker(T,A{1},A{2},A{3}))
%!test % 1d
%! T = randn(2,1);
%! A = randn(3,2);
%! out = tucker(T,A);
%! ref = A*T;
%! assert(out,ref,1e-13)
%!test % 2d
%! T = randn(2,3);
%! A{1} = randn(3,2);
%! A{2} = randn(4,3);
%! out = tucker(T,A);
%! ref = A{1}*T*A{2}.';
%! assert(out,ref,1e-13)
%!test % 3d
%! T = randn(2,3,4);
%! A{1} = randn(3,2);
%! A{2} = randn(4,3);
%! A{3} = randn(5,4);
%! out = tucker(T,A);
%! ref = mump(mump(mump(T,A{1},1),A{2},2),A{3},3);
%! assert(out,ref,1e-13)
%!test % 4d
%! T = randn(2,3,4,5);
%! A{1} = randn(3,2);
%! A{2} = randn(4,3);
%! A{3} = randn(5,4);
%! A{4} = randn(6,5);
%! out = tucker(T,A);
%! ref = mump(mump(mump(mump(T,A{1},1),A{2},2),A{3},3),A{4},4);
%! assert(out,ref,1e-13)
%!test % tensorization
%! A{1} = randn(2,1);
%! A{2} = randn(3,1);
%! A{3} = randn(4,1);
%! A{4} = randn(5,1);
%! out = tucker(1,A);
%! ref = tensorize(A);
%! assert(out,ref,1e-13)
%!test %tensor with implicit last dimension
%! T = randn(2,3,4);
%! A{1} = randn(2);
%! A{2} = randn(3);
%! A{3} = randn(4);
%! A{4} = randn(5,1);
%! out = tucker(T,A);
%! ref = mump(mump(mump(mump(T,A{1},1),A{2},2),A{3},3),A{4},4);
%! assert(out,ref,1e-13)
%!test % complex
%! T = randn(2,3,4)+1i*randn(2,3,4);
%! A{1} = randn(3,2)+1i*randn(3,2);
%! A{2} = randn(4,3)+1i*randn(4,3);
%! A{3} = randn(5,4)+1i*randn(5,4);
%! out = tucker(T,A);
%! ref = mump(mump(mump(T,A{1},1),A{2},2),A{3},3);
%! assert(out,ref,1e-13)
%!test % Jump some modes
%! T = randn(2,3,4,5);
%! A1 = randn(3,2);
%! A2 = randn(4,3);
%! A3 = randn(5,4);
%! A4 = randn(6,5);
%! out = tucker(T,[],A2,A3,A4);
%! ref = mump(mump(mump(T,A2,2),A3,3),A4,4);
%! assert(out,ref,1e-13)
%! out = tucker(T,A1,[],A3,A4);
%! ref = mump(mump(mump(T,A1,1),A3,3),A4,4);
%! assert(out,ref,1e-13)
%! out = tucker(T,A1,A2,[],A4);
%! ref = mump(mump(mump(T,A1,1),A2,2),A4,4);
%! assert(out,ref,1e-13)
%! out = tucker(T,A1,A2,A3,[]);
%! ref = mump(mump(mump(T,A1,1),A2,2),A3,3);
%! assert(out,ref,1e-13)
%! out = tucker(T,[],A2,A3,[]);
%! ref = mump(mump(T,A2,2),A3,3);
%! assert(out,ref,1e-13)
%! out = tucker(T,A1,[],[],A4);
%! ref = mump(mump(T,A1,1),A4,4);
%! assert(out,ref,1e-13)
%! out = tucker(T,A1,[],A3,[]);
%! ref = mump(mump(T,A1,1),A3,3);
%! assert(out,ref,1e-13)
%! out = tucker(T,[],[],A3,[]);
%! ref = mump(T,A3,3);
%! assert(out,ref,1e-13)
%!error
%! tucker();
%!error
%! tucker(randn(3));
%!error
%! tucker(randn(3),[]);
