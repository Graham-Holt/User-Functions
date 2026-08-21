function Ainv = invFF(A,modulo)
% invFF(A,modulo) computes the finite field matrix inverse
% 
% Graham Holt, August 2026. Updated August 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% invFF(A,modulo)
% Ainv = invFF(A,modulo)
% 
%% Description
% invFF(A,modulo) returns a square matrix B such that B*A = A*B = I in the
% finite field "modulo"

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