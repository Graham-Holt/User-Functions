function P = intpartition(n,m,leq)
% intpartition(n) generates the integer partitions of "n".
% 
% Graham Holt, June 2026. Updated June 2026.
% Embry-Riddle Aeronautical University
% 
%% Syntax
% intpartition(n)
% P = intpartition(n)
% 
%% Description
% intpartition(n) returns a cell array with the integer partitions of "n".
%
% intpartition(n,m) returns a cell array with the integer partitions of "n"
% containing at most "m" summands.
%
% intpartition(n,m,1) returns a cell array with the integer partitions of
% "n" containing exactly "m" summands.

p = zeros(1,n); p(1) = n;
P = {}; k = 1;
if nargin < 3
    leq = 0;
end

while true
    if ((k < m)&&(~leq)) || (k == m)
       P{end+1} = p(1:k); 
    end
    if k==n
        return;
    end

    r = k - find(p>1,1,'last') + 1;
    k = find(p>1,1,'last');
    p(k) = p(k) - 1;

    while r>p(k)
        p(k+1) = p(k);
        r = r - p(k);
        k = k + 1;
    end

    p(k+1) = r;
    k = k + 1;
end