function E = pagemexp(M)
% pagemexp(M) performs matrix exponential across array pages
% 
% Graham Holt, March 2026. Updated May 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% pagemexp(M)
% L = pagemexp(M)
% 
%% Description
% pagemexp(M) takes in a 3D array of matrices and performs the matrix
% exponential across each page of the input

E = zeros(size(M));
for k = 1:size(M,3)
    E(:,:,k) = expm(M(:,:,k));
end