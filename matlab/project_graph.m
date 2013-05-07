function varargout = project_graph(v, A, AA, L, D, P)
% PROJECT_GRAPH    Project a point onto the graph of a linear operator.
%
%   If A is an m x n matrix and v = [v1;v2] is a vector such that v1 is of
%   length n and v2 is of length m, then 
%
%     z = project_graph(v,A);
%   
%   gives the projection of v onto the graph { (x,y) | y = Ax } of A. The
%   vector z can then be sliced into its first n components and the remaining
%   m components.
%
%   If A is dense and m <= n, then the projection involves computing
%   the Cholesky factor L of I + A*A'. One can obtain A*A' and L by calling
%
%     [z AA L] = project_graph(v,A)
%
%   and can pass in AA and L to avoid them being recomputed via
%
%     z = project_graph(v,A,AA,L)
%
%   This allows for factorization caching when calling project_graph
%   multiple times for the same A.
%
%   This applies similarly when m >= n, only AA is A'*A and L is the
%   Cholesky factor of I + A'*A.
%
%   This also applies similarly when A is sparse, in which case we
%   use a permuted LDL factorization of some matrix. In this case,
%   one can obtain P, L, and D via
%
%     [z L D P] = project_graph(v,A);
%
%   and can pass these back in via
%
%     z = project_graph(v,A,[],L,D,P);
%
%   Note that L has different meanings depending on whether or not A is sparse.

    [m n] = size(A);
    c = v(1:n);
    d = v(n+1:end);

    if issparse(A)
        if ~exist('P','var') || isempty(P) || ...
                ~exist('L','var') || isempty(L) || ~exist('D','var') || isempty(D)
            K = [ speye(n) A' ; A -speye(m) ];
            [L,D,P] = ldl(K);
        end

        z = P * (L' \ (D \ (L \ (P' * sparse([ c + A'*d ; zeros(m,1) ])))));

        varargout(1) = {z};
        varargout(2) = {L};
        varargout(3) = {D};
        varargout(4) = {P};
    else
        if m <= n
            if ~exist('AA','var') || isempty(AA)
                AA = A*A';
            end
            if ~exist('L','var') || isempty(L)
                L = chol(eye(m) + AA);
            end
            y = L \ (L' \ (A*c + AA*d));
            x = c + A'*(d - y);
        else
            if ~exist('AA','var') || isempty(AA)
                AA = A'*A;
            end
            if ~exist('L','var')
                L = chol(eye(n) + AA);
            end
            x = L \ (L' \ (c + A'*d));
            y = A*x;
        end

        varargout(1) = {[x;y]};
        varargout(2) = {AA};
        varargout(3) = {L};
    end

end
