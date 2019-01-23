%% info:

%% 

% (t/c)_r = s the thickness to chord ratio at the planform root at the
% centerline of the aircraf;


% thickness to chord ratio at the planform root at the centerline of the aircraft
t_c_at_root = 8.7/100 ;

% is the ratio of the thickness to chord ratios at the root and the tip of the planform

% NOTE: This airplane has uniform cross section through out the whole thing
% so the ratio will be 1!
Tau = t_c_at_root / t_c_at_root ;

% is the taper ratio which is the ratio of the root chord to the tip chord of the planform


% the followin numbers are estimation based on visual inspection!

Ct = 0.05;  % cord length @ tip, estimated via visual inspection and functions match

Cr = 0.23;  % Cord length @ root, you can see it from pictures
Lambda = Ct/Cr ; 


Cfe = 0.004; %  equivelent skin friction drag


% is the exposed planform area of the planform not including the area within the fuselage

Sexp_plf = 0.63 ;

Swet_wing = 2 * Sexp_plf * ( 1 + 0.25 * (t_c_at_root) * ( ( 1 + (Tau*Lambda)) / ( 1 + Lambda)) ) * 3.7 ;

