% Final Glider Design
clc; 
clear all; 
% Equation to go from m to feer and to inches 

meter = .2462;
inches = meter * 3.28084 * 12 

% Overall Glider configuration variables

% chord length of wing at planform root
Cr =  .0889;

% Chrod length at the tip
Ct = Cr/2;

% Thickness of wing at planform root
% .052 is the % thickness of the airfoil
Percent_thickness = .121;
T_r =  Cr*(Percent_thickness);


% Thickness of wing at tip
T_t = Ct*(Percent_thickness);


% wing span
b = .6; % m



%% Wing Area

% Trapezoid with offset



S_ref = b*(Cr+Ct)/2 ;

% Ellipse 
%S_ref = (Cr*b*pi) / 2 ;


% Aspect Ratio
AR = b^2 / S_ref;


% Sweep Angle
SW = atan((Cr-Ct)/b);

% Taper Ratio (Lamda)
TR = Ct/Cr ;

% tau is the ratio of thickness to chord ratiors at the root and the tip of
% the planform 
tau = (T_r/Cr)/(T_t/Ct) ; 



%% Calculate the wetted area 


% Calculate the wetted area of the wings 

S_esp_plf = S_ref*2 ; 


S_wet_plf = 2 * S_esp_plf * (1 + (.25*(T_r/Cr))* (1+(tau*TR))/(1+TR));


% Total wetted area

%S_wet = .125; 
S_wet = .032 + 2*(.043) + 2*(.008) + .014;

%% Weight Calculations

% define density of materials

rho_Balsa = 160; %kg/m^3
rho_Bamboo = 780;
rho_DepronFoam = 40; %kg/m^3
rho_Coroplast = .9 * .001* (1/(1*10^-6)); %kg/cm^3
rho_PlyWood = 600; %kg/m^3
rho_PosterBoard = 250*.001; %kg/m^3
rho_Tissue = 24*.001; %kg/m^30
rho_ABS = 1250;
rho_Pinkfoam = 24.8;

% fuselage weight
%V_fue = (.014 + .009)*.0254;
V_fue = 2.422*10^-4;

W_fuselage = V_fue * rho_Pinkfoam* 9.81;

% Ballast weight

W_ballast = .001 * 9.81;


% tail weight
% 0.0016 m is the thickness of the balsa wood

V_HT = 1.031*10^-5;

W_tail = 2 * V_HT *rho_Balsa* 9.81 ;


% Using the Selig S1223 high lift low Reynolds number airfoil
% Estimating the area of the airfoil to get an approximation of the Volume
% of the wing

%V_wing = S_ref*T_t;Area_airfoil_tip*Thickness_CrossSection 

%V_wing =  1.97*10^-5; 


Area_airfoil_root = 8.386*10^-5;
Area_airfoil_tip = Area_airfoil_root/2;
Area_airfoil_middle = Area_airfoil_root * .75;
% area of the airfoil at the root multiplied by the thickness of each cross section
Thickness_CrossSection = .0015875; % m (1/16 inches)

V_structure = (Area_airfoil_root)*Thickness_CrossSection + Area_airfoil_tip*Thickness_CrossSection + Area_airfoil_middle*Thickness_CrossSection + Area_airfoil_root*.625*Thickness_CrossSection + Area_airfoil_root*.875*Thickness_CrossSection;

%W_wing = 2*(V_structure * rho_Balsa + ((2*S_ref*(1+(.25*(T_r/Cr)*((1+tau+TR)/(1+TR))))) * rho_Tissue))*9.81;

W_wing = 0.0147*2;
% New equation for the Weight of the wing

%W_wing = V_airfoil*rho_Balsa


%% Calculate the Center of Gravity 

%cg_desired = .01;
cg_desired = 0.2*Cr;

% input all the distances of the center of gravity of each component to the
% 
BALL =  .069;

FUSE = .099-.069;

VSTAB = .318-.069;

HSTAB = .333-.069;

% Ballast

x_cgBallast = BALL + cg_desired;

% Wing 

x_cgWing = .4*Cr - cg_desired;

% Fuselage 

x_cgFuselage = FUSE - cg_desired;

% Vertical Stabilizer is included in the fuselage
x_cgVstabilizer = VSTAB - cg_desired;

% Horizontal Stabilizer 

x_cgHstabilizer = HSTAB - cg_desired;


% center of gravity for the whole plane 
%M_cg = (x_cgBallast * W_ballast) + (x_cgWing * W_wing) + (x_cgFuselage * W_fuselage) + (x_cgHstabilizer * W_tail);

%M_cg = (x_cgWing * W_wing) + (x_cgFuselage * W_fuselage) + (x_cgHstabilizer * W_tail) ;

M_cg = (x_cgWing * W_wing) + (x_cgFuselage * W_fuselage) + (x_cgHstabilizer * W_tail) ;
% Calculate the necessary ballast weight to get the cg at 20% of the chord
DesiredBallastWeight = abs((M_cg/x_cgBallast));


M_cgNew = (x_cgWing * W_wing) + (x_cgFuselage * W_fuselage) + (x_cgHstabilizer * W_tail) + (DesiredBallastWeight*x_cgBallast);
% Calculate the total weight with the desired ballast weight 
Wtape = (.02*9.81);
Wdowelrod = (.01*9.81);

TotalWeight = (DesiredBallastWeight + W_fuselage + W_tail + W_wing + Wtape + Wdowelrod );


%TotalMass = (W_ballast + W_fuselage + W_tail + W_wing)/9.81;
TotalMass = (TotalWeight)/9.81 ;

Wing_loading = TotalWeight / S_ref;


%% Calculate Cl and required wing loading 


max_range = 100; % meters
launch_height = 9.5 ; % meters 

Cl_Cd_maxRange = max_range/launch_height;

Cl_Cd_max = 20.847;

% Reynolds number
[T, a1, P, rho_altitude] = atmosisa(1632);

Velocity = 5; % m/s
length_fuselage = .338 % m
% mu at 5000 ft (from data package)
mu_altitude = 1.74*(10^-5); % 
% rho_atmosphere slug/ft^3


Re = (rho_altitude * Velocity * length_fuselage) / mu_altitude;

% skin friction coefficient 
S_fe = .074/Re^0.2;

% zero lift drag 
Cd_0 = S_fe*(S_wet/S_ref);

%  

Cl = Cl_Cd_max * 2 * Cd_0; 

% Cl_Cd_actual = Cl/(2*Cd_0); 

% Calculate the wing loading required based on Cl value and Velocity
% required (3-7 m/s the goal is 3 max 7)

%Velocity_required = sqrt(Wing_loading / (Cl * .5 * rho_altitude*S_ref));

%% Plot the wing loading vs velocity

V_inf = linspace(0,7,500);

WingLoading = Cl*.5*rho_altitude*(V_inf.^2);
Wing_Loading = ( .062*9.81) / S_ref ;

plot(V_inf,WingLoading)
hold on
yline(Wing_loading);
hold off
xlabel('Velocity (m/s)');
ylabel('Wing Loading');
title(' Velocity vs Wing Loading');
legend('V vs Wing Loading','Calculated Wing Loading');


%% Calculate the Longitudinal Static Stability and Trim and Lateral Directional Stability

% Longitudinal Static Stability
Chord = linspace(Ct,Cr,1000);
ChordAvg = sum(Chord)/1000;

S_h = 0.003*2; % meters squared 


VH = (S_h * x_cgHstabilizer) / (S_ref*ChordAvg);

VH_flight = (S_h*.26) / ((S_ref)*ChordAvg);

% Lateral Directional Stability 
S_v = 0.004;% meters squared 

VV = (S_v * x_cgVstabilizer) / (S_ref*b);

DesiredValue = .035/(S_ref*b);


%% Calculate the angle of attack

a_0 = (1.868-1.3)/5 ; 

% given value of e
e = .9 ;

k = 1/(pi*e*AR);

% angle of attack where lift = 0 
AoA_l0 = -11 ;

a1 = a_0 / (1 + ((57.3*a_0)/(pi*e*AR)));

AoA = (Cl/a1) - AoA_l0;


Cd = Cd_0 + (k * Cl^2);


%Cl_Cd_maxActual = Cl/Cd;

%% Dihedral Angle

B = 5;

r = (b * Cl * B)/ (x_cgVstabilizer);

%% Cl Horizontal Tail

CL_HT = (.055 - Cl*( (cg_desired/Cr) - (Cr*.25/Cr) ) )/-VH ;


%% Masses
fus = .01;
SingleWing = 0.003;
SingleStabilizer = 8.981*10^-4;
span_wing = .8;


%% Stuff for Flight Day

% PreFlight_weight = ;
% Ground_range = ;
% t_flight = ;
% 
% 
% Flight_LD_max = ;
% Alpha 

height = 6.5;
a = 0.1074;
alphaL_0 = -2.75; 


S_flight = zeros(4,1);
V_trim = zeros(4,1);
CL = zeros(4,1);
alpha = zeros(4,1);

Ground_range = [7, 6.9, 7.8, 16];




PreFlight_weight = [.063*9.81, .0631*9.81, .0634*9.81, .064 *9.81];


t_flight = [2,3,2,4];





for i = 1:4
    
S_flight(i) = sqrt(height^2 + Ground_range(i)^2);

V_trim(i) = S_flight(i)/t_flight(i) ; 

CL(i) = PreFlight_weight(i)/ (.5 * rho_altitude * V_trim(i)^2 * S_ref);

alpha(i) = (CL(i) / (a)) + alphaL_0 ;

end



V_trimActual = mean(V_trim);

L_D_flight = Ground_range / height;

plot(alpha,L_D_flight,'*')

L_D_flight = L_D_flight';
f = fit(alpha,L_D_flight,'poly3','Normalize','on','Robust','Bisquare');
plot(f,alpha,L_D_flight)
title('L/D vs alpha Experimental');
xlabel('alpha (degrees)');
ylabel('L/D');







