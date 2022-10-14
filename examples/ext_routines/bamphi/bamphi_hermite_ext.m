function x = bamphi_hermite_ext( x, z, conjugated, max_No )
% HERMITESORT extend the interpolation set x with points from z so that these
% are sorted "a' la Hermite", which is an improved version of the Leja sorting.
%
% Input:
%      - x, given interpolation set;
%      - z, set where we pick the extensions;
%      - conjugated, is true if we desire to pick the points from z in complex
%        conjugated pairs, when possible;
%      - max_No, is the maximum number of points that we pick from z. Useful for
%        performing approximation to Leja sequences over peculiar bounded regions.
% Output:
%      - x, extended set.
%
% REFERENCE: M.Caliari, P.Kandolf and F.Zivcovich, Backward error analysis of
% polynomial approximations for computing the action of the matrix exponential,
% Bit Numer. Math. 58 (2018), 907, https://doi.org/10.1007/s10543-018-0718-9 .
%
%      Author: Franco Zivcovich
%       email: franco.zivcovich@gmail.com
% Last Update: 18 August 2020.
%
  % Check inputs
  if nargin < 2
    error('Not enough input arguments');
  end
  if ( nargin < 3 ) || isempty( conjugated )
    conjugated = false;
  end
  if isempty( x )
    [ ~, idx ] = max( abs( z ) );
    x( 1 ) = z( idx( 1 ) ); z( idx( 1 ) ) = [];
    if ( conjugated && imag( x ) )
      idx = find( z == conj( x ) );
      if ~isempty( idx )
        x( 2 ) = z( idx( 1 ) ); z( idx( 1 ) ) = [];
      end
    end
    clear idx;
  end
  % Start body
  n_x = length( x ); n_z = length( z );
  if ( nargin < 4 ) || isempty( max_No )
    max_No = n_x + n_z;
  end
  % Create w_x = [ 1; -sum(out); ...; prod(out); zeros( k,1 )]
  w_x = eye( 1 + max_No, 1 );
  for j = 1 : n_x
    w_x( 2 : j+1 ) = w_x( 2 : j+1 ) - x( j ) * w_x( 1 : j );
  end
  % Set identification: build matrix K of knots and matrix M of multiplicities
  K = zeros( n_z, 1 ); % over-allocate K
  M = false( n_z, max_No ); % over-allocate M
  k = 0; id_max = 0;
  while ~isempty( z ),
    % search and destroy
    k = k + 1;
    K( k ) = z( 1 );
    id = ( z == K( k ) );
    z( id ) = [];
    % count
    s_id_x = sum( x == K( k ) );
    s_id_z = sum( id );
    M( k, 1 + s_id_x : s_id_x + s_id_z ) = true;
    id_max = max( id_max, s_id_x + s_id_z );
  end
  % Sorting core
  for i = 1 : id_max
    [ x, w_x ] = hermite_sort( x, w_x, i - 1, K( M( 1:k,i ) ), conjugated, max_No );
  end
end
function [ x, w_x ] = hermite_sort( x, w_x, dif, z, conjugated, max_No )
  n_x = length( x ); n_z = length( z ); N = n_x + n_z;
  % Build rectangular Vandermonde
  % X = [ cumprod( repmat( z, 1, N - dif - 1 ), 2, 'reverse' ), ones( n_z,1 ) ];
  X = [ fliplr( cumprod( repmat( z, 1, N - dif - 1 )' )' ), ones( n_z,1 ) ];
  % Differentiation vector ( dif_vec * w_x = d_difw_x )
  dif_vec = exp( gammaln( N : -1 : dif+1 ) - gammaln( N - dif : -1 : 1 ) ).';
  while ( ~isempty( z ) && n_x <= max_No )
    n_x = n_x + 1; n_z = n_z - 1;
    % find
    [ ~, idx ] = max( abs( X( :,n_z+1:N-dif ) * ( dif_vec( 1:n_x-dif ) .* w_x( 1:n_x-dif ) ) ) );
    % set and destroy
    x( n_x ) = z( idx );
    z( idx ) = [];
    X( idx,: ) = [];
    if conjugated
      idx = ( z == conj( x( n_x ) ) );
      if any( idx )
        n_x = n_x + 1; n_z = n_z - 1;
        % set and destroy
        x( n_x ) = z( idx );
        z( idx ) = [];
        X( idx,: ) = [];
        % refresh w_x
        w_x( 2:n_x+1 ) = w_x( 2:n_x+1 ) - 2 * real( x( n_x ) ) * w_x( 1:n_x ) ...
                         + [ 0; abs( x( n_x ) )^2 * w_x( 1:n_x-1 ) ];
        continue,
      end
    end
    % refresh w_x
    w_x( 2:n_x+1 ) = w_x( 2:n_x+1 ) - x( n_x ) * w_x( 1:n_x );
  end
end
%endfunction
