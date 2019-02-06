clear
clc


%% Estimation of the surface area of the UAV (Exluding wings)

% Assume cylinder extending from front fuselage to 0.45 m
h = .45;
r = 0.08;
A1 = pi * r * (r + sqrt(h^2 + r^2));

% Assume cylinder extending from 0.45 m to front of horizontal stabilizer
h2 =0.9800;
A2 = 2 * pi * r * h;

% Verical stabilizer 1 - Assume triangle
base = .28;
height = .4;
A3 = 2 * .5 * base * height;

% Verical stabilizer 2 - Assume rectangle
l = 0.13;
b = .4;
A4 = 2 * l * b;

% Estimate the horizontal stabilizer
Length = 0.75;
Width = 0.13;

A5= 2 * Length * Width;

% Sum of areas
ATotal  = A1 + A2 + A3 + A4 +A5
