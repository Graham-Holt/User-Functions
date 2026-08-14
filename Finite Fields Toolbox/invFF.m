function Ainv = invFF(A,modulo)

if size(A,1)~=size(A,2)
    error('Matrix must be square.');
end
if mod(det(A),modulo)==0
    error('Matrix is singular.');
end

n = size(A,1);
if n==1

    [g,Ainv,~] = gcd(A,modulo);
    
    if g~=1
        Ainv = [];
    else
        Ainv = mod(Ainv,modulo);
    end
    
else
    
    Ainv = rrefFF([A eye(n)],modulo);
    Ainv = Ainv(:,(n+1):end);

end