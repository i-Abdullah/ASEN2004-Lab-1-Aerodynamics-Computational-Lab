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
Tau = (8.7/100)*0.23 / ((30/100)*0.23) ;

% is the taper ratio which is the ratio of the root chord to the tip chord of the planform


% the followin numbers are estimation based on visual inspection!

Ct = 0.363636 * 0.23 ;  % cord length @ tip, estimated via visual inspection and functions match

Cr = 0.23;  % Cord length @ root, you can see it from pictures
Lambda = Ct/Cr ;




Cfe = 0.004; %  equivelent skin friction drag


% is the exposed planform area of the planform not including the area within the fuselage

Sexp_plf = 0.63 ;

Swet_wing = 2 * Sexp_plf * ( 1 + 0.25 * (t_c_at_root) * ( ( 1 + (Tau*Lambda)) / ( 1 + Lambda)) ) * 3.7 ;



%% Estimated surface area of airplane form front fuselage to tip of wing
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

%%

Atotal = AFuselage+Swet_wing

Sref = 0.63 ; %surface refrence area


CdMin = (Atotal/Sref) * Cfe
