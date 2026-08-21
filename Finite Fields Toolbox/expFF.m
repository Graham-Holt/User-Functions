function n = expFF(b,a,modulo)
% expFF(b,a,modulo) computes exponentiation in finite fields
% 
% Graham Holt, August 2026. Updated August 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% expFF(b,a,modulo)
% n = expFF(b,a,modulo)
% 
%% Description
% expFF(b,a,modulo) performs exponentiation in the finite field "modulo"

b = int64(b);
a = int64(a);

if a<0

    a = -a;
    for k = 1:length(b(:))
        b(k) = invFF(b(k),modulo);
    end

end
if ~isinteger(a)

    warning('Rational approximation of exponent used.');

    [aN,aD] = rat(a);
    b = expFF(b,int64(aN),modulo);

    n = expFF(0:(modulo-1),int64(aD),modulo);
    I = find(n==b);

    if isempty(I)
        error('Answer does not exist.');
    elseif ~isscalar(I)
        warning('Answer is not unique');
    end

    n = n(I(1));

else

    n = ones(size(b));
    for k = 1:a
        n = mod(b.*n,modulo);
    end

end