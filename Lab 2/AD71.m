%% Info: this script is done for the airplane AD71, A power-off glider

% part of CU's Lab, SP 18, 2004: Aircraft design.

% note:

% x = center of mass for that element
%% Define constants.

clear
clc
close all

%% material info

% using 1/8'th of inch balsa wood

thickness_balsa = (1/32) * 0.0254 ;  
rhoBalsa = 160 ; % kg/m^3
g = 9.81;
%%

BoulderAlt = 1624; % in m
LaunchAlt = 9 ; % 
Viscosity = 1.74e-5 ;


% get info using atmosisa function
[ Temp SpeedSound Pressure Density ] = atmosisa(BoulderAlt+LaunchAlt);

% Wing span

Span = 0.7 ; % in meters.
TaperRatio = 0.5; %
ChordRoot = 0.05;
ChordTip = ChordRoot*TaperRatio ;

MAC = (( ChordTip + ChordRoot ) / 2 ); % Mean aerodynamic chord.

WingArea = Span * MAC ;

AspectRatio = ((Span)^2)/WingArea;

Velocity = 6;

ReynoldsWing = ( Density * Velocity * MAC ) / Viscosity ; % Reynold numbers for wings

SweepAngle = 0; % creates drag, more lift, we don't need it.

WingSurfaceArea = 0.195 ; % m^2;
WingSurfaceArea = 0.1323 ; % m^2;


WingVolume = WingSurfaceArea * thickness_balsa ;
WingMass = WingVolume*rhoBalsa ;

WingWeight = WingMass * g ;


%% Stabilizers:

% HS = Horizontal stablizer.
% VS = Vertical stablizer.

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


% VS Is just a trapzoidal!
VS_span = Span/6; % span of vs is 3rd 
VS_root = 0.06 ; % length at root.
VS_tip = 0.06*0.5; % tip length.


VS_Area = ( ( VS_root + VS_tip ) / 2 ) * VS_span ; 

% thickness * surface area to  ge volume
VS_Volume = VS_Area * thickness_balsa ; 

VS_Mass = VS_Volume*rhoBalsa ;
VS_Weight = VS_Mass*g ;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

HS_root = VS_tip*4; % attached to fueslage, same like VS height in
HS_tip = 0.05; %area at outside edge or length at outside edge
HS_span = Span/6 ; %from side of fueslage to the outside, length to the outside
HS_Area = ( ( HS_root + HS_tip ) / 2 ) * HS_span ; 

% thickness * surface area to  ge volume
HS_Volume = HS_Area * thickness_balsa ; 

HS_Mass = HS_Volume*rhoBalsa ; % in kg.
HS_Weight = HS_Mass*g ;
HS_Weight = 2*HS_Weight;


%% fueslage

Fuselage_length = 0.2; % across the body
Fueslage_width = 0.003; % up or down

Fueslage_Area = Fuselage_length*Fueslage_width ;
Fueslage_Volume = Fueslage_Area*thickness_balsa ;
Fueslage_Mass = Fueslage_Volume*rhoBalsa ;
Fueslage_Weight = Fueslage_Mass*g;

xcg_Fueslage = Fuselage_length/2 ; 

%% cg

%centroid of trapazoid in x direction, let c be symmetric on both sides

% http://www.efunda.com/math/areas/trapezoid.cfm

c_VS = sqrt(((VS_root-VS_tip)^2 + (VS_span^2))) ;
xcg_VS = ( (2*VS_tip*c_VS) + (VS_tip)^2 + (c_VS*VS_root) + (VS_root)^2 ) / ( 3 * ( VS_root+VS_tip)) ;


c_HS = sqrt(((HS_root-HS_tip)^2 + (HS_span^2))) ;
xcg_HS = ( (2*HS_tip*c_HS) + (HS_tip)^2 + (c_HS*HS_root) + (HS_root)^2 ) / ( 3 * ( HS_root+HS_tip)) ;

% location of thee wing: half way on the fueslage

xcg_Fueslage = Fuselage_length/2 ;

% location of thee wing: half way on the fueslage
x_wing_leadingEdge_location = xcg_Fueslage;

xcg_wanted = (27/100)*ChordRoot  ; 

% HS Starts at VS_root which is right at the end of fueslage.
xcg_HS = (Fuselage_length - VS_root) + xcg_HS - xcg_wanted - x_wing_leadingEdge_location ;
xcg_VS = (Fuselage_length - VS_root) + xcg_VS - xcg_wanted - x_wing_leadingEdge_location ;
xcg_Fueslage = (Fuselage_length/2) - xcg_wanted - x_wing_leadingEdge_location ;
xcg_wing = ((40/100) * ChordRoot) + x_wing_leadingEdge_location - xcg_wanted - x_wing_leadingEdge_location ;




%% total weight and mass, and moment

% total weight
Weight_tot = VS_Weight + HS_Weight + Fueslage_Weight + WingWeight ;
% total mass
Mass_tot = VS_Mass + HS_Mass + Fueslage_Mass + WingMass ;

%total area
Sref = 2*HS_Area + VS_Area + WingSurfaceArea + Fueslage_Area ;


%total moment
TotalMoment = (xcg_VS*VS_Weight) + (xcg_HS*HS_Weight) + (xcg_wing*WingWeight) + (xcg_Fueslage*Fueslage_Weight) ;

%IMPORTANT: Relative to leading edge of root!
xcg_aircraft = ( TotalMoment / Weight_tot ) + x_wing_leadingEdge_location ; 

% this moment will be at distance of length of fueslage from begining up
% until the root of wing, which is trnaslated into mass of:

MassNeeded = TotalMoment / (x_wing_leadingEdge_location+xcg_wanted) ; 

fprintf('Mass Needed to balance the airplane %0.3f', MassNeeded);
fprintf(' kg \n');

fprintf('total mass %0.3f', MassNeeded+Mass_tot);
fprintf(' range required (0.04 to 0.08) kg \n');

%% Horizontal tail volume coefficient:

VH = ( 2*HS_Area * ( xcg_HS - xcg_wanted ) ) / ( Sref*MAC) ;

fprintf('Horizontal tail volume coefficient %0.3f ', VH);
fprintf(' ideal range (0.3 to 0.6) \n');

%% wing loading

WingLoading = Weight_tot/Sref ;