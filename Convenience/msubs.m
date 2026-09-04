function funcM = msubs(varargin)

func = varargin{1};
if isscalar(symvar(func)) && nargin < 3
    var = symvar(func);
    M = varargin{2};
else
    var = varargin{2};
    M = varargin{3};
end

n = size(M,1);
[V,J] = jordan(M);

L = diag(J);
[L,~,I] = unique(L);

funcJ = [];
for k = 1:length(L)
    nBlock = sum(I==k);
    funcJk = zeros(nBlock);
    for j = 1:nBlock
        dFuncL = subs(diff(func,j-1),var,L(k));
        funcJk = funcJk + diag(dFuncL*ones(nBlock-j+1,1)/factorial(j-1),j-1);
    end
    funcJ = blkdiag(funcJ,funcJk);
end
funcM = V*funcJ/V;