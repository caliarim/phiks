function W = sylvphi(A,U,p,APhi1,AAPhi2)
% Phi1 or phi2 functions using Sylvester method
% INPUT:
% A: cell of two matrices
% U: matrix
% p: either 1 or 2. If 1 it computes phi1, if 2 it computes phi2
% OUTPUT:
% if p=1 unvec of phi1(I \otimes A{1} + A{2} \otimes A{2})U(:) 
% if p=2 unvec of phi2(I \otimes A{1} + A{2} \otimes A{2})U(:) 
%
% Code extracted from https://github.com/jmunoz022/Kronecker_EI
% and it works just in 2D

  if (nargin <= 3)
      APhi1{1} = A{1}*phipade(A{1},1);
      APhi1{2} = A{2}*phipade(A{2},1);
  end
  if(p==1)
    W = APhi1{1}*(U*APhi1{2}.') + U*APhi1{2}.' + APhi1{1}*U;
    W = sylvester(A{1},A{2}.',W);
  elseif(p==2)
    if (nargin <= 4)
      AAPhi2{1} = A{1}*(A{1}*phipade(A{1},2));
      AAPhi2{2} = A{2}*(A{2}*phipade(A{2},2));
    end
    W = APhi1{1}*(U*APhi1{2}.');
    W = W+AAPhi2{1}*U+U*AAPhi2{2}.';
    W = sylvester(A{1},A{2}.',W);
    W = sylvester(A{1},A{2}.',W);
  else
    error('The code computes either phi1 (p=1) or phi2 (p=2)')
  end
end
%!test
%! x = linspace(0,1,4)';
%! y = linspace(1,2,5)';
%! [X,Y] = ndgrid(x,y);
%! U = X.*Y;
%! A{1} = toeplitz([-2,1,0,0]);
%! A{2} = toeplitz([-2,1,0,0,0]);
%! K = kron(speye(5),A{1})+kron(A{2},speye(4));
%! ref1 = phi1m(K)*U(:);
%! ref2 = phi2m(K)*U(:);
%! assert(sylvphi(A,U,1)(:),ref1,1e-10)
%! assert(sylvphi(A,U,2)(:),ref2,1e-10)
%! APhi1{1} = A{1}*phi1m(A{1});
%! APhi1{2} = A{2}*phi1m(A{2});
%! AAPhi2{1} = A{1}*(A{1}*phi2m(A{1}));
%! AAPhi2{2} = A{2}*(A{2}*phi2m(A{2}));
%! assert(sylvphi(A,U,1)(:),ref1,1e-10)
%! assert(sylvphi(A,U,2)(:),ref2,1e-10)
%! assert(sylvphi(A,U,1,APhi1)(:),ref1,1e-10)
%! assert(sylvphi(A,U,2,APhi1,AAPhi2)(:),ref2,1e-10)
