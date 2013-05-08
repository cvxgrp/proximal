function x = project_affine_set(v, A, b)
% PROJECT_AFFINE_SET    Project a point into an affine set.
%
%   project_affine_set(v,A,b) is the projection of v onto
%   the affine set { x | Ax = b }.

    x = v - pinv(A)*(A*v - b);
end
