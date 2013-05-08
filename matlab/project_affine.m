function [x C d] = project_affine(v, A, b, C, d)
% PROJECT_AFFINE    Project a point into an affine set.
%
%   project_affine(v,A,b) is the projection of v onto
%   the affine set { x | Ax = b }.
%
%   You can also call the function as
%
%     [x C d] = project_affine(v,A,b);
%
%   and then call it again with different argument v2 via
%
%     x2 = project_affine(v2,A,b,C,d);
%
%   If calling the function repeatedly with the same A and b,
%   all evaluations after the initial one will be much faster
%   when passing in the cached values C and d.

    if exist('C','var') && ~isempty(C) && exist('d','var') && ~isempty(d)
        x = C*v + d;
    else
        % x = v - pinv(A)*(A*v - b);
        pA = pinv(A);
        C = eye(length(v)) - pA*A;
        d = pA*b;
        x = C*v + d;
    end
end
