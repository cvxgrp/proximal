function x = prox_cvx(v, lambda, f, l, u)
% PROX_CVX    The proximal operator of a generic function.
%
%   prox_cvx(v,lambda,f) is the proximal operator of f with parameter lambda
%   evaluated at v. Here, f is a closure, possibly described by an anonymous
%   function. It is also possible to provide lower and upper bounds, so
%   prox_cvx(v,f,l,u) is the prox operator of f + I_[l,u].
%
%   For example,
%
%     prox_cvx(v, lambda, (@(x) norm(x,1))
%
%   is equivalent to the soft thresholding operator.
%
%   WARNING: This is a *very* inefficient way of evaluating a proximal
%   operator. This function is mainly for rapid prototyping and
%   testing custom implementations of particular proximal operators.
%   It may also be useful when suffering from extreme laziness.


    if ~exist('l', 'var') || isnan(l)
        l = -Inf;
    end

    if ~exist('u', 'var') || isnan(u)
        u = Inf;
    end

    L = 1/(2*lambda);
    [m n] = size(v);

    if min(m,n) == 1
        cvx_begin quiet
            variable x(max(m,n))
            minimize(f(x) + L*sum_square(x - v))
            subject to
                v >= l;
                v <= u;
        cvx_end
    elseif m == n && m > 1
        cvx_begin quiet
            variable x(n,n) symmetric
            minimize(f(x) + L*pow_pos(norm(x - v, 'fro'), 2))
            subject to
                x == semidefinite(n)
        cvx_end
    else
        x = nan(size(v));
    end
end
