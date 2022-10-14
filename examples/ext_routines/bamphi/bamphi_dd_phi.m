function [ dd, F, s ] = bamphi_dd_phi( x, l, L )
% dd_phi_ compute phi_l’s divided differences at x.
%
% Input:
%      -  x, interpolation points;
%      -  l, index of the phi_l function desired;
%      -  L, number of consecutive points equal to zero in the tail of z (for
%            performance only, default L = 0).
% Output:
%      - dd, divided differences of phi_l at x.
%
% REFERENCE: Zivcovich, F. Fast and accurate computation of divided differences
% for analytic functions, with an application to the exponential function.
% Dolomites Research Notes on Approximation, 12(1), 28-42, (2019).
%
%      Author: Franco Zivcovich
%       email: franco.zivcovich@gmail.com
% Last Update: 18 August 2020.
%
  if ( ( nargin < 3 ) || isempty( L ) ), L = 0; end
  if ( ( nargin < 2 ) || isempty( l ) ), l = 0; end
  n = length( x ) + l;
  x = [ zeros( 1,l ), reshape( x, [ 1, n-l ] ) ];
  % Compute scaling
  F = repmat( x, n, 1 ); F = tril( F - F.' );
  s = max( ceil( max( abs( F( : ) ) ) / 3.5 ), 1 );
  % Compute F_0
  N = n + 30;
  dd = [ 1, 1 ./ cumprod( ( 1 : N ) * s ) ];
  for j = n - L : -1 : 1
    dd( j+L : n-1 ) = dd( j+L : n-1 ) + x( n-L : -1 : j+1 ) .* dd( j+L+1 : n );
    for k = N : -1 : n
      dd( k ) = dd( k ) + x( j ) * dd( k + 1 );
    end
  end
  F( n-L+1 : n, n-L+1 : n ) = toeplitz( dd( 1 : max( L,1 ) ) ); % max( L,1 ) because toeplitz([]) triggers an error
  for j = n - L : -1 : 1
    for k = n - j : -1 : 2
      dd( k ) = dd( k ) + F( k+j, j ) * dd( k+1 );
    end
    F( j,j+1 : n ) = dd( 2 : n-j+1 );
  end
  F( 1 : n+1 : n^2 ) = [ exp( x( 1 : n-L ) / s ), ones( 1,L ) ];
  F = triu( F );
  % Square F_0 into F_s
  dd = F( 1,: );
  for k = 1 : s - 1
    dd = dd * F;
  end
  dd = dd( l+1 : n ).';
end
