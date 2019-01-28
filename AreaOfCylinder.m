clear 
clc

% Estimated surface area of airplane form front fuselage to tip of wing
r = 0.0800; % Radius in m
h = 0.45; % Hight/length of cone in m
A = pi * r * ( r + sqrt( h^2 + r^2)); % Area from front fuselage to tip of wing in m^2
BottomArea = pi * r^2; % Compute bottom area of the cone in m^2
A1 = A - BottomArea; % Substract bottom are of the cone in m^2

% Estimated surface area of airplane from front wing to to front canard
h2 = 1.56 - h - 0.13; % Height/length of cone in m 
A2 = pi * r * ( r + sqrt( h2^2 + r^2)); % Area of bottom of cone in m^2
A3 = A2 - BottomArea; % Substract bottom of wing from total area in m^2

% Total surface area of the airplane excluding wings
AFuselage = A3 + A1 % Total area of fuselage

Atotal = AFuselage+4.7634








