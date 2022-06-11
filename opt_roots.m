function r = opt_roots(k)
  function [w,f,g,ngp] = vars(r, x)
    p = prod(1 - x' ./ r, 1)';
    w = x ./ (1 - p.^2);
    f = sqrt(w) .* p;
    q = sum(1 ./ (x' - r), 1)';
    g = x .* (1./(2*w) + q);
    ngp = sum((p'.^2 + r./(x'-r))./(x'-r), 1)';
  end

  r = .5 - .5*cos([1:k]'/(k+.5) * pi);
  x = .5 - .5*cos((.5+[1:k-1]')/(k+.5) * pi);

  dr = r; drsize = 1;
  while drsize > 128*eps
    dx = x; dxsize = 1;
    while dxsize > 128*eps
      dxsize = norm(dx, inf);
      [~,~,g,ngp] = vars(r,x);
      dx = g ./ ngp; x = x + dx;
    end
    x1 = [x;1]; [w,f,~,~] = vars(r,x1);
    f0 = sqrt(.5/sum(1./r));
    J = f0^3./r'.^2 + w.*abs(f)./(r'.*(x1-r'));
    drsize = norm(dr, inf);
    dr = -J \ (f0 - abs(f)); r = r + dr;
  end
end
