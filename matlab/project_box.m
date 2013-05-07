function x = project_box(v, l, u)
% PROJECT_BOX    Project a point onto a box (hyper-rectangle).
%
%   project_box(v,l,u) is the projection of v onto
%   the set { x | l <= x <= u }.

    x = max(l, min(v, u));
end
