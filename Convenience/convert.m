function x2 = convert(x1,u1,u2,abs)
% convert(x1,u1,u2) converts the input from the first set of units to the
% second
% 
% Graham Holt, Feb. 2024. Updated May 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% convert(x1,u1,u2)
% convert(___,abs)
% x2 = convert(__)
% 
%% Description
% convert(x1,u1,u2) takes in the units "u1" and "u2" as strings and 
% converts "x1" to "x2" using those units.
% 
% convert(x1,u1,u2,abs) uses the binary indicator "abs" to determine if
% temperature units represent absolute temperature or differential
% temperature
% 
% List of Units: https://www.mathworks.com/help/symbolic/units-list.html

u1 = str2symunit(u1); u2 = str2symunit(u2);

if abs
    x2 = double(separateUnits(unitConvert(x1*u1,u2,'Temperature','absolute')));
else
    x2 = double(separateUnits(unitConvert(x1*u1,u2,'Temperature','difference')));
end