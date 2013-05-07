function x = project_exp(v)
% PROJECT_EXP    Project points onto the exponential cone.
%
%   When v is a 1 x 3 vector, project_exp(v) is the projection
%   of v onto the exponential cone. When v is an n x 3 vector,
%   project_exp(v) projects each row of v onto the cone
%   in a vectorized fashion.
%
%   For reference, the exponential cone and its dual are given by
%     Kexp   = { (x,y,z) | ye^(x/y) <= z, y > 0 }
%     Kexp^* = { (u,v,w) | u < 0, -ue^(v/u) <= ew } cup { (0,v,w) | v,w >= 0 }

    r = v(:,1); s = v(:,2); t = v(:,3);
    x = nan(size(v));

    % v in cl(Kexp)
    idx = ( (s.*exp(r./s) <= t & s > 0) | (r <= 0 & s == 0 & t >= 0) );
    x(idx,:) = v(idx,:); 

    % -v in Kexp^*
    idx = ( (-r < 0 & r.*exp(s./r) <= -exp(1).*t) | (-r == 0 & -s >= 0 & -t >= 0) );
    x(idx,:) = 0;

    % special case with analytical solution
    idx = (r < 0 & s < 0);
    x(idx,:) = v(idx,:);
    x(idx,2) = max(x(idx,2),0);
    x(idx,3) = max(x(idx,3),0);

    % minimize ||x - v||^2 subject to se^{r/s} = t via primal-dual Newton method
    % these components are computed serially, so much slower
    idx = find(isnan(x(:,1)));

    g     = @(w)   w(2)*exp(w(1)/w(2)) - w(3);
    gradg = @(w)   [ exp(w(1)/w(2)); exp(w(1)/w(2))*(1 - w(1)/w(2)); -1 ];

    alpha = 0.001; beta = 0.5;

    for i = 1:length(idx)
        disp('newton')
        u = v(idx(i),:)';
        u(2) = max(u(2),1); u(3) = max(u(3),1);
        y = 1; % dual variable

        r = @(w,z) [ w - v(idx(i),:)' + z*gradg(w); g(w) ];

        for iter = 1:100
            KKT = [ eye(3)+y*hessg(u), gradg(u) ; gradg(u)', 0 ];
            z = KKT \ -r(u,y);
            du = z(1:3);
            dy = z(4);

            % backtracking line search
            t = 1;
            ustep = u + t*du; ystep = y + t*dy;
            while ustep(2) < 0 || (norm(r(ustep, ystep)) > (1 - alpha*t)*norm(r(u, y)))
                t = beta*t;
                ustep = u + t*du; ystep = y + t*dy;
            end

            u = ustep; y = ystep;

            if abs(g(u)) < 1e-8 && norm(r(u,y)) <= 1e-8
                x(idx(i),:) = u;
                break;
            end
        end
    end
end

function h = hessg(w)
    r = w(1); s = w(2); t = w(2);
    h = exp(r/s)*[ 1/s,    -r/s^2,   0;
                   -r/s^2, r^2/s^3,  0;
                   0,      0,        0 ];
end
