function P = intpartition(n)
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

p = zeros(1,n); p(1) = n;
P = {}; k = 1;

while true
    P{end+1} = p(1:k);
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