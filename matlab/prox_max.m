function x = prox_max(v, lambda)
% PROX_MAX    The proximal operator of the max function.
%
%   prox_max(v,lambda) is the proximal operator of the max
%   of the entries of v. This function is not vectorized.


    TOL      = 1e-8;
    MAX_ITER = 100;
    rho = 1/lambda;

    n = length(v);
    tl = min(v) - 1/n;
    tu = max(v);

    g = @(t) (sum(max(0, rho*(v - t))) - 1); 

    iter = 0;
    while tu - tl > TOL && iter < MAX_ITER 
        t0 = (tl + tu)/2;
        if sign(g(t0)) == sign(g(tl))
            tl = t0; 
        else
            tu = t0; 
        end 
        iter = iter + 1;
    end 

    x = min(t0,v);

    if tu - tl > TOL
        warning('Algorithm did not converge.');
    end

end
