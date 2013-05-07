function x = project_pos(v)
% PROJECT_POS    Project a point onto the nonnegative orthant.
%
%   project_pos(v) is the positive part of v.

    x = max(v,0);
end
