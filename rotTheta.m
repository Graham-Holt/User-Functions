function R = rotTheta(theta,u)
% rotTheta(theta,u) generates a rotation matrix according to the Principle 
% Axis Theorem
% 
% Graham Holt, March 2026. Updated May 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% rotTheta(theta,u)
% R = rotTheta(theta,u)
% 
%% Description
% rotTheta(theta,u) takes in an angle "theta" (in radians) and a 3D unit
% vector "u" and returns the matrix which rotates a 3D vector about "u" by
% "theta".

if length(u)~=3
    error('Input must be columnm vector with length 3');
end

R = expm(theta*crossmat(u/norm(u)));