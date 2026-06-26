function T = rootTree(p,s)

if p==1
    T = {ones(s+1,1)};
    return;
end

A = sym('a%d%d',[s s]); A = blkdiag(tril(A,-1),1);
K = intpartition(p-1); k = 1;
for j = 1:length(K)
    Phi = {}; M = [];
    for i = 1:length(K{j})
        Phi{i} = rootTree(K{j}(i),s);
        M(i) = size(Phi{i},2);
    end
    I = ones(1,length(Phi));
    while all(I<=M)
        T{k} = [ones(s,1); 1/p];
        for m = 1:length(Phi)
            T{k} = T{k}.*(A*Phi{m}(:,I(m)));
        end

        k = k + 1;
        if all(I==M)
            break;
        end

        I(end) = I(end)+1;
        while any(I>M)
            I(find(I>M,1,'last')+1) = I(find(I>M,1,'last')+1)+1;
            I(find(I>M,1,'last')) = 1;
        end
    end
end

function parts = intpartition(n,s)
%
% DISCLAIMER: This function was licensed via the MathWorks File Exchange.
% Credit to the original writer of the script given here:
% https://www.mathworks.com/matlabcentral/fileexchange/33616-integer-partitions?s_tid=ta_fx_results
%
%INTPARTITIONS performs integer partition, i.e. the partition of of a set 
%containing homogenous elements. The function generates a cell array 
%containing a list of vectors representing all possible ways 
%of partitioning a set containing intIn number of identical elements 
%without order awareness. The numerical representation in the 
%output describes the partitions as: {[3 1]} = [1 1 1 | 1].
%   
%   Syntax:
%   intpartition(n)
%   intpartition(n,s)
%
%   The resulting partions can also be seen as writing the input n as a sum
%   of positive integers. An optional argument s can be supplied to output
%   a subset of partitions with number of parts less than or equal to s.
%
%   Number of ways of partitioning is according to sequence:
%   http://oeis.org/A000041 - (1), 1, 2, 3, 5, 7, 11, 15, 22, 30, 42, ...
%
%   Example 1: intpartition(4) gives {[1 1 1 1],[1 1 2],[1 3],[2 2],4}
%   Example 2: intpartition(10,2) gives {[3,7],[4,6],[5,5],10}
%

n = round(n);
if ~isscalar(n)
    error('Invalid input. Input must be a scalar integer')
end
if ~exist('s','var')
    s = inf;
end
    
parts={};
a = zeros(1,n);
k = 1;
y = n - 1;
while k ~= 0
    x = a(k) + 1;
    k = k-1;
    while 2*x <= y
        a(k+1) = x;
        y = y - x;
        k = k + 1;
    end
    l = k + 1;
    while x <= y
        a(k+1) = x;
        a(l+1) = y;
        
        if s==inf || k+2==s
            parts(end+1) = {a(1:(k+2))};
        end
        x = x + 1;
        y = y - 1;
    end
    a(k+1) = x + y;
    y = x + y - 1;
    if s==inf || k+1==s
        parts(end+1) = {a(1:(k+1))};
    end
end