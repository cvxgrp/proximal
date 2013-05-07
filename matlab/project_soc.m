function x = project_soc(v)
% PROJECT_SOC    Project a point onto the second-order cone.
%
%   Suppose v is an (n+1)-dimensional vector, so v = [t;v0], where t is a
%   scalar and v0 is n-dimensional. Then project_soc(v) is the projection of
%   (t,v0) onto the second-order cone; the result is also (n+1)-dimensional.

% assume v = (t, v0)
    nv = norm(v(2:end));
    x  = nan(size(v));
    if nv <= -v(1)
        x = zeros(size(v));
    elseif nv <= v(1)
        x = v;
    else
        r = 0.5*(1 + v(1)/nv);
        x(1)     = r*nv;
        x(2:end) = r*v(2:end);
    end
end
