%% Compute the minimum of a given function using Newton-Raphson iteration.
%
% fun     : Scalar function to be minimized
% grad    : Gradient function of 'fun'
% hessian : Hessian function of 'fun'
% x0      : Initial vector
% tol     : Tolerance
% maxiter : Maximum iteration
function x = newton(fun, grad, hessian, x0, tol, maxiter)
  x = x0;
  eps = 1.0;
  iter = 0;
  while (eps > tol && iter < maxiter)
    y = fun(x);
    dfdx = grad(x);
    H = hessian(x);
    dx = -H\dfdx;
    x = x + dx;
    eps = norm(dx);
    iter = iter + 1;
  end
end