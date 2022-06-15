function beta = opt_coef(k)
  x = cos([1:k]'*pi/(k+.5));  % Quadrature nodes
  w = (1-x)/(k+.5);  % Quadrature weights

  % 4th kind Chebyshev polynomials W evaluated at x
  W = zeros(k, k);
  W(:, 1) = 1;
  if k >= 2; W(:, 2) = 2*x+1; end
  for i=3:k
    W(:,i) = (2*x) .* W(:,i-1) - W(:,i-2);
  end

  r = opt_roots(k);  % Roots of optimal polynomial p
  lambda = (1-x)/2;  % Nodes transformed to [0, 1]
  p = prod(1-lambda'./r, 1)';  % p(lambda)
  alpha = W' * (w .* p);
  beta = 1 - cumsum((2*[0:k-1]'+1) .* alpha, 1);
end
