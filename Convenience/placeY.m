function K = placeY(A,B,C,D,eig)

[n,m] = size(B);
[p,~] = size(C);
s = length(eig) - n;
if s < 0
    error('Too few eigenvalues to place.');
elseif s > n-m-p+1
    s = n-m-p+1;
    eig = eig(1:(n+s));
end

A = blkdiag(A,zeros(s));
B = blkdiag(B,eye(s));
C = blkdiag(C,eye(s));
D = blkdiag(D,zeros(s));

K = sym('k%d%d',[m p]+s);

S = solve(charpoly(A - B*K/(eye(p+s) + D*K)*C) == poly(eig), symvar(K));

for k = 1:(m+s)
    for j = 1:(p+s)
        K(j,k) = S.("k"+num2str(k)+num2str(j));
    end
end
K = double(K);