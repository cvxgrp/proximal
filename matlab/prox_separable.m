function [x iter] = prox_separable(v, lambda, fp, l, u, x0, tol, MAX_ITER)
% PROX_SEPARABLE   Evaluate the prox operator of a fully separable function.
%
% Arguments:
%
%  v is the point at which to evaluate the operator.
%  fp is a subgradient oracle for the function.
%  lambda (optional) is the proximal parameter; defaults to 1.
%  l (optional) is a lower bound for x; defaults to -Inf.
%  u (optional) is an upper bound for x; defaults to Inf.
%  x0 (optional) is a value at which to warm start the algorithm.
%  tol (optional) is a stopping tolerance.
%  MAX_ITER (optional) is the maximum number of iterations.
%
% Examples:
%
%  v = randn(n,1);
%  x = prox_separable(v, 1, @(w) sign(w));
%  [x iter] = prox_separable(v, 1, @(w) sign(w));
%
% This function can be called in vectorized form if fp is vectorized,
% i.e., if fp works elementwise, then v, l, u, and x0 can be vectors.

    n = length(v);

    if ~exist('lambda', 'var') || isnan(lambda) || isempty(lambda)
        lambda = 1;
    end
    rho = 1/lambda;

    if ~exist('l', 'var') || any(isnan(l)) || isempty(l)
        l = -inf(n,1);
    end
    if ~exist('u', 'var') || any(isnan(u)) || isempty(u)
        u = inf(n,1);
    end
    if ~exist('x0', 'var') || any(isnan(x0)) || isempty(x0)
        x0 = zeros(n,1);
    end
    if ~exist('tol', 'var') || isnan(tol) || isempty(tol)
        tol = 1e-8;
    end
    if ~exist('MAX_ITER', 'var') || isnan(MAX_ITER) || isempty(MAX_ITER)
        MAX_ITER = 500;
    end

    iter = 0;
    x = max(l, min(x0, u));

    while any(u-l > tol) && iter < MAX_ITER
        g = fp(x) + rho*(x - v);

        idx = (g > 0); 
        l(idx) = max(l(idx), x(idx) - g(idx)/rho);
        u(idx) = x(idx);

        idx = ~idx;
        u(idx) = min(u(idx), x(idx) - g(idx)/rho);
        l(idx) = x(idx);

        x = (l + u)/2;
        iter = iter + 1;
    end

    if any(u-l > tol)
        fprintf(2, 'Warning: %d entries did not converge; max interval size = %f.\n', ...
            sum(u-l > tol), max(u-l));
    end
end
