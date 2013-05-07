function x = prox_quad(v, lambda, A, b)
% PROX_QUAD    The proximal operator of a quadratic.
%
%   prox_quad(v,lambda,A,b) 

    rho = 1/lambda;
    n = size(A);
    if issparse(A)
        x = (A + rho*speye(m)) \ (rho*v - b);
    else
        x = (A + rho*eye(m)) \ (rho*v - b);
    end
end
