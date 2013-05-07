function x = prox_sum_square(v, lambda)
% PROX_SUM_SQUARE    Proximal operator of sum-of-squares.
%
%   prox_sum_square(v,lambda) is the proximal operator of
%   (1/2)||.||_2^2 with parameter lambda.

    x = (1/(1 + lambda))*v;
end
