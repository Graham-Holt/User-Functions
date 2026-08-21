function a = logFF(b,n,modulo)
% logFF(b,n,modulo) computes logarithms in finite fields
% 
% Graham Holt, August 2026. Updated August 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% logFF(b,n,modulo)
% a = logFF(b,n,modulo)
% 
%% Description
% logFF(b,n,modulo) performs logarithms in the finite field "modulo"

for k = 0:(modulo-2)
    B = expFF(b,k,modulo);

    if B == n
        a = k;
        return;
    end

    if k > 0 && B == 1
        error('Logarithm does not exist.');
    end
end
error('Logarithm does not exist.');