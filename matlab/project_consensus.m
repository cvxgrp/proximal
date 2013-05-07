function x = project_consensus(varargin)
% PROJECT_CONSENSUS    Bring a set of points into consensus.
%
%   project_consensus(v1,v2,...,vn) returns the elementwise
%   average of v1, ..., vn, i.e., the projection of the vi
%   onto the consensus set.

    N = length(varargin);
    x = zeros(size(varargin{1}));
    for i = 1:N
        x = x + varargin{i};
    end
    x = x./N;
end
