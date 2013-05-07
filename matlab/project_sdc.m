function X = project_sdc(V)
% PROJECT_SDC    Project a point onto the semidefinite cone.
%
%   project_sdc(V), where V is a symmetric n x n matrix, is the
%   projection of V onto the semidefinite cone.

    [W,D] = eig(V);
    X = W*max(0,D)*W';
end
