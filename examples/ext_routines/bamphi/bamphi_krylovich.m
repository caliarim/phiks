function varargout = bamphi_krylovich( str, jump, Atilde, A, f, p, id, rly, opts, info )

  persistent H V n_V

  switch str

    case 'build'

      low_stor = ( opts.low_stor && opts.iom );

      % Initialize data
      V{ 1 } = f;
      V( opts.r_arn + p + 1 : length( V ) ) = [];
      V_ = flip( eye( p,1 ) );

      % Preliminaries
      for k = 1 : p
        idk    = low_stor + k * not( low_stor );
        idk_p1 = low_stor + ( k + 1 ) * not( low_stor );
        [ V{ idk_p1 }, V_( :, idk_p1 ) ] = Atilde( V{ idk }, V_( :, idk ) );
      end
      idk_p1 = low_stor + ( p + 1 ) * not( low_stor );
      n_V = sqrt( V{ idk_p1 }' * V{ idk_p1 } ); %n_V
      V{ idk_p1 } = V{ idk_p1 } / n_V;

      % Perform ( iom ) - Arnoldi
      H = zeros( opts.r_arn + 1 );
      for j = 1 : opts.r_arn

        idpj_m1 = ( mod( j - 2, 3 ) + 1 ) * low_stor + ( p + j - 1 ) * not( low_stor );
        idpj    = ( mod( j - 1, 3 ) + 1 ) * low_stor + ( p + j     ) * not( low_stor );
        idpj_p1 = ( mod( j + 0, 3 ) + 1 ) * low_stor + ( p + j + 1 ) * not( low_stor );

        V{ idpj_p1 } = A( V{ idpj } );
        for k = j : -1 : max( ( j - 1 ) * ( abs( info.A.skew ) || opts.iom ),1 )

          idk = ( mod( k - 1, 3 ) + 1 ) * low_stor + ( p + k ) * not( low_stor );

          H( k,j ) = V{ idk }' * V{ idpj_p1 };
          V{ idpj_p1 } = V{ idpj_p1 } - H( k,j ) * V{ idk };
        end
        H( j+1,j ) = sqrt( V{ idpj_p1 }' * V{ idpj_p1 } );
        if ( H( j+1,j ) == 0 ) % check for convergence
          break,
        end
        if info.A.skew
          H( j,j+1 ) = info.A.skew * H( j+1,j );
        end
        V{ idpj_p1 } = V{ idpj_p1 } / H( j+1,j );
      end

      % Return
      varargout{ 1 } = H;

    case 'assemble'

      ret = 1;
      m = length( V ) - 1;
      d = info.A.his.d{ id }{ jump };
      x = info.A.his.x{ id };

      % Determine coefficients for Arnoldi space
      [ cf_out, cf_aux ] = cof_krylov( info.A.his.ts( id ), H, x, d, p, n_V );
      if rly
        cf_out = real( cf_out );
      end

      w_ = cf_out( p : -1 : 1 ); % = V_( :, 1:p ) * cf_out( 1:p );
      e_t = min( ceil( 1.1 * m ), opts.m_r_arn );
      w = 0; w1 = 0;
      for j = 1 : length( V )
        w  = w + cf_out( j ) * V{ j };
        if j > p
          w1 = w1 + cf_aux( j - p ) * V{ j };
        end
        if ( abs( cf_out( j ) ) < opts.tol * sqrt( w' * w ) )
          e_t = min( e_t, j );
        end
      end
      % Wrap up interpolation at Ritz's values
      w1 = ( info.A.his.ts( id ) * d( m + 2 ) ) * ( A( w1 ) - x( m + 1 ) * w1 );
      w = w + w1;
      % Now take care of the extended interpolation points
      if rly
        x = real( x );
        d = real( d );
        w1 = real( w1 );
        w = real( w ); w_ = real( w_ );
      end
      c1 = norm( w1,inf );
      for j = 3 : length( x ) - m
        w1 = ( info.A.his.ts( id ) * d( m + j ) / d( m + j - 1 ) ) * ( A( w1 ) - x( m + j - 1 ) * w1 );
        w = w + w1;
        c2 = norm( w1,inf );
        % [ c1, c2, norm( w, Inf ) * opts.tol ], pause
        if ( ( c1 + c2 ) < norm( w, Inf ) * opts.tol )
          ret = 0;
          break;
        end
        c1 = c2;
      end
      its = m + j - 1;

      % Return
      varargout{ 1 } = w;
      varargout{ 2 } = w_;
      varargout{ 3 } = its;
      varargout{ 4 } = e_t;
      varargout{ 5 } = ret;

  end

end

%
function [ f_out, f_aux ] = cof_krylov( t, H, x, d, p, n_V )
% Compute in a competitive time the vector f_out = expm( t*H ) * e1;

  % Compute phipm( H ) * eye( length( H ), 1 )
  m = length( H ) - 1;
  f_out = [ d( p + 1 ); zeros( m,1 ) ];
  f_aux = [ t * ( H( 1 ) - x( p + 1 ) ); t * H( 2 ); zeros( m-1,1 ) ];
  for j = 2 : m
    f_out = f_out + f_aux * d( p + j );
    f_aux = t * ( H * f_aux - x( p + j ) * f_aux );
  end
  f_out = f_out + f_aux * d( p + m + 1 );

  % Now adjust coefficients
  f_out = [ [ 1; cumprod( t * ones( p - 1, 1 ) ) ] .* d( 1 : p ); ( n_V * t^p ) * f_out ];
  f_aux = ( n_V * t^p ) * f_aux;

end
%endfunction
