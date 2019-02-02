%% info:

%% 

clear
clc


% (t/c)_r = s the thickness to chord ratio at the planform root at the
% centerline of the aircraf;


% thickness to chord ratio at the planform root at the centerline of the aircraft
t_c_at_root = 8.7/100 ;

% is the ratio of the thickness to chord ratios at the root and the tip of the planform

% NOTE: This airplane has uniform cross section through out the whole thing
% so the ratio will be 1!
Tau = (8.7/100) / (8.7/100) ;

% is the taper ratio which is the ratio of the root chord to the tip chord of the planform


% the followin numbers are estimation based on visual inspection!

Ct = 0.363636 * 0.23 ;  % cord length @ tip, estimated via visual inspection and functions match

Cr = 0.23;  % Cord length @ root, you can see it from pictures
Lambda = Ct/Cr ;




Cfe = 0.004; %  equivelent skin friction drag


% is the exposed planform area of the planform not including the area within the fuselage

Sexp_plf = 0.63 ; % only wing platform area
Non_exposed = 0.23*0.16 ; % the non exposed side of the wing that attaches to the fueslage.

Total_Exposed = Sexp_plf - Non_exposed;

Swet_wing = 2 * Total_Exposed * ( 1 + 0.25 * (t_c_at_root) * ( ( 1 + (Tau*Lambda)) / ( 1 + Lambda)) )  ;



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

%%
Atotal = ATotal + Swet_wing
Sref = 0.63 ; %surface refrence area


CdMin = (Atotal/Sref) * Cfe
