function GammaK = fov_contour(extreigs,npoints)
%
% function GammaK = fov_contour(extreigs,npoints)
%
% extreigs = [SR,LR,SI,LI]
if (nargin == 1)
  npoints = 50;
end
if (extreigs(2) == extreigs(1)) % vertical line
  if (extreigs(4) == extreigs(3)) % a single point
    GammaK = extreigs(1)+1i*extreigs(3);
  else
    GammaK = extreigs(1)+1i*linspace(extreigs(3),extreigs(4),npoints);
  end
elseif (extreigs(4) == extreigs(3)) % horizontal line
  GammaK = linspace(extreigs(1),extreigs(2),npoints)+1i*extreigs(3);
else % rectangle
  hmx = max(ceil(npoints/2*(extreigs(2)-extreigs(1))/...
                 (extreigs(2)-extreigs(1)+extreigs(4)-extreigs(3))),2);
  hmy = max(ceil(npoints/2*(extreigs(4)-extreigs(3))/...
                 (extreigs(2)-extreigs(1)+extreigs(4)-extreigs(3))),2);
  GammaK = linspace(extreigs(1),extreigs(2),hmx)+1i*extreigs(3);
  GammaK = [GammaK(2:end-1),extreigs(2)+1i*linspace(extreigs(3),extreigs(4),hmy)];
  GammaK = [GammaK(1:end-1),linspace(extreigs(2),extreigs(1),hmx)+1i*extreigs(4)];
  GammaK = [GammaK(1:end-1),extreigs(1)+1i*linspace(extreigs(4),extreigs(3),hmy)];
end
%!demo
%! GammaK = fov_contour([-2,2,-2,2]);
%! plot(real(GammaK),imag(GammaK),'*')
%! length(GammaK)
%! axis equal
%!demo
%! GammaK = fov_contour([-100,0,-1,1]);
%! plot(real(GammaK),imag(GammaK),'*')
%! length(GammaK)
%! axis equal
%!demo
%! GammaK = fov_contour([-2,2,-2,2],100);
%! plot(real(GammaK),imag(GammaK),'*')
%! length(GammaK)
%! axis equal
%!demo
%! GammaK = fov_contour([-2,2,2,2]);
%! plot(real(GammaK),imag(GammaK),'*')
%! length(GammaK)
%! axis equal
%!demo
%! GammaK = fov_contour([2,2,-2,2]);
%! length(GammaK)
%! plot(real(GammaK),imag(GammaK),'*')
%! axis equal
%!demo
%! GammaK = fov_contour([0,0,0,0]);
%! length(GammaK)
%! plot(real(GammaK),imag(GammaK),'*')
%! axis equal
