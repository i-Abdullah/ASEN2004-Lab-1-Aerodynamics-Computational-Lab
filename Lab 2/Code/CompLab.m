%% info


%% Housekeeping

clear;
clc;
close all;

%% Read data

%CFD Data
TempsetCFD = xlsread("Tempest UAS CFD flight data CFD.xlsx");
AlphaCFD = TempsetCFD(:,1);
CL_CFD = TempsetCFD(:,2);
CD_CFD = TempsetCFD(:,3);


Airfoil2D_Data = xlsread("xf-e71-il-50000-n5.xls");
% Airfoil2D_Data = xlsread("xf-s1223-il-50000-n5.xls");
% Remove uneeded cells
Airfoil2D_Data(1:8,:) = [];
Alpha2D = Airfoil2D_Data(:,1);
Cl_2D = Airfoil2D_Data(:,2);
Cd_2D = Airfoil2D_Data(:,3);

%% Find angle and Cm at L = 0
MinCl = min(abs(Cl_2D));
for i = 1:length(Cl_2D)
    if Cl_2D(i) == -MinCl
       AOA_Zero_Lift = Alpha2D(i);
       Cm = Airfoil2D_Data(i,5);
    end  
end
WingSpan = .6;  % m
S        = .04; % m^2
AspectRatio = WingSpan^2/S; 
EfficRatio = 0.9;

%% find where lift is == 0 

% for this we will take the avg between the negative and possitive

% first positive number
Postiv = find(Cl_2D>0,1);

% Alpha (AOA) when lift is == 0 
%AOA_Zero_Lift = mean([Alpha2D(Postiv);Alpha2D(Postiv-1)]);
%Alpha0 = (Cl_2D(Postiv) - Cl_2D(Postiv - 1))/(Alpha2D(Postiv) - Alpha2D(Postiv - 1)) * (Alpha2D(Postiv));
% pick point in the middle

[ r c ] = size(Cl_2D);
MidPoint = ceil(r/2);


% The slope of the linear region for the 2d wing.
a0 = (Cl_2D(MidPoint+1)-Cl_2D(MidPoint))/((Alpha2D(MidPoint+1)-Alpha2D(MidPoint)));




% the lift curve slope for 3D
a3D_Wing_liftCurveSlope = (a0)/(1+ ( ( 57.3 * a0 ) / ( pi * EfficRatio * AspectRatio) )) ;

CL_3DWing_Estimated = a3D_Wing_liftCurveSlope .* ( Alpha2D - AOA_Zero_Lift);

%% Induced drag, and  wing drag

InducedDrag_Wing3D = (CL_3DWing_Estimated).^2 ./ (pi*EfficRatio*AspectRatio);

% the wing drag 
WingDrag = Cd_2D + InducedDrag_Wing3D;

% get alpha when wing drage is min to be used in equation 3.5;

Alpha_wing_mindD = Alpha2D(find(min(WingDrag)==WingDrag));

%% the whole aircraft drag

% find where the the minumm drage happens on the whole wing, get CL
% corresponding. 

% The equation used is 3.4a in the pdf, it's on page 4
% the component of the equation is in the next page.

e0 = 1.78 * ( 1 - 0.045 .* AspectRatio.^(0.68)) - 0.64 ;
k1 = 1 / ( pi * e0 * AspectRatio ) ;

%CL_whenDragIsMin_for_the_whole_airplane
CL_MinD_Airplane = a3D_Wing_liftCurveSlope * ( Alpha_wing_mindD - AOA_Zero_Lift);
CD_Min_theWhole_Airplane = 0.012599210884700 ;


%%%%%%%%%%%%%%%%%%%%% Start of AM's changes %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Surface area of all components - not including ballast 
S_Wing       = 2 * .0432;  % m^2
S_Fuselage   = .032;        % m^2
S_VS         = .014;        % m^2
S_HS         = 2 * .008;    % m^2
Sref         = .01 * 4;     % m^2

% Wetted areas 
Swet_Tot     = S_Wing + S_Fuselage + S_VS + S_HS; % m^2
Swet_Tot     = 0.148;

% Coefficient of friction 
[T, a, P, rho] = atmosisa(1633);
V              = 7;                           % m/s
L_Fuselage     = .337;                          % m
u              = 1.74*10^(-5);                % kg/m s
Re             = (rho * V * L_Fuselage) / u;  % Unitless
Cfe            = 0.074 / Re^(.2);             % Unitless
 
% Drag of whole airplane at minimum lift
CD_Min_theWhole_Airplane = Cfe * Swet_Tot / Sref % Unitless
%%%%%%%%%%%%%%%%%%%%%%%%%% End of AM's changes %%%%%%%%%%%%%%%%%%%%%%%%%


ParasiteDrag_CD0_wholeAirplane = CD_Min_theWhole_Airplane + k1*(CL_MinD_Airplane)^2 ;

k2 = -2*k1*CL_MinD_Airplane;

CD_theWholeAirplane_Polar = ParasiteDrag_CD0_wholeAirplane + k1*(CL_3DWing_Estimated).^2 + k2*(CL_3DWing_Estimated);



%% Calculating the dihedral angle
% B = 5;
% Xacv = 
% r = (WingSpan * CL * B) / (Xacv - Xcg); 


%% L/D : is it 3d just wing or the whoel aircraft/


L_D_WholeAirplane = CL_3DWing_Estimated ./ CD_theWholeAirplane_Polar ; 
%L_D_WholeAirplane = 16 / 6.4;

L_D_CFD = (CL_CFD./CD_CFD);

%% velocity to achieve max range and max endurance

GOTA = 6.4; % Kg, groos weight
GOTAWeight = GOTA*9.81;
Density = 1.0324 ; %kg/m^3 @ 1.8 km.
WingArea = 0.63 ; % wing area.
V_MaxRangeEndurance_Equation = @(CL_V) sqrt ( (2 *( GOTAWeight/WingArea)) / ((Density)*CL_V));

CL_Max_Range = sqrt( ParasiteDrag_CD0_wholeAirplane/k1);
CL_Max_Endurance = sqrt( (3*ParasiteDrag_CD0_wholeAirplane)/k1);


CL_CD_Ratio_Max_Endurance = GOTAWeight/ParasiteDrag_CD0_wholeAirplane;
CL_CD_Ratio_Max_Range = GOTAWeight/ParasiteDrag_CD0_wholeAirplane;

V_Max_Range = V_MaxRangeEndurance_Equation(CL_Max_Range);
V_Max_Endurance = V_MaxRangeEndurance_Equation(CL_Max_Endurance);


%based on the assumption that CL of the wing is CL for the airplane, we can estimate the AOA from graph of equations

Vstall = sqrt((2 * 9.81 * 0.063)/(Sref * 1.0442 * max(CL_3DWing_Estimated)))


AOA_Max_Range = (CL_Max_Range/a3D_Wing_liftCurveSlope)+AOA_Zero_Lift;
AOA_Max_Endurance = (CL_Max_Endurance/a3D_Wing_liftCurveSlope)+AOA_Zero_Lift;
max(L_D_WholeAirplane)
%% plot

figure(1)

plot(Alpha2D,CL_3DWing_Estimated,'*-','LineWidth',1)
hold on
% plot(AlphaCFD,CL_CFD,'*-','LineWidth',1)
% hold on
plot(Alpha2D,Cl_2D,'*-','LineWidth',1)
hold on
refline(0)
hold off

% legend('3D Calculated','CFD','2D','Location','SouthEast')
legend('3D Calculated','2D','Location','SouthEast')
xlabel('\alpha \circ ')
ylabel('Coefficient of Lift')
title('Lift Curve Comparison')
grid minor


% - - -


figure(2)

plot(Cl_2D,Cd_2D,'*-','LineWidth',1)
hold on
% plot(CL_3DWing_Estimated,WingDrag,'*-','LineWidth',1)
% hold on
refline(0)
hold off

legend('finite wing drag','Location','NorthWest')
% legend('infinite wing drag','finite wing drag','Location','NorthWest')
xlabel(' Coefficients of Lift ')
ylabel(' Coefficient of Drag')
title('Drag Polar Comparison: The wings')
grid minor

% - - - 

figure(3)

plot(Alpha2D,L_D_WholeAirplane,'*-','LineWidth',1)
% hold on
% plot(AlphaCFD,L_D_CFD,'*-','LineWidth',1)
hold on
plot(9.3265,1.0769,'*')
hold on
plot(24.8862,1.0615,'*')
hold on
plot(8.0073,1.2,'*')
hold on
plot(12.2638,2.4615,'*')
refline(0)
hold off

legend('3D Wings full L/D','Location','NorthWest')
% legend('3D Wings full L/D','CFD L/D','Location','NorthWest')
xlabel(' \alpha \circ ')
ylabel('L/D')
title('L/D Comparison for Glider Design')
grid minor

% - - - 

figure(4)

plot(CL_3DWing_Estimated,CD_theWholeAirplane_Polar,'*-','LineWidth',1)
% hold on
% plot(CL_CFD,CD_CFD,'*-','LineWidth',1)
hold on
refline(0)
hold off

legend('estimated Polar Drag','Location','NorthWest')
% legend('estimated Polar Drag','estimated CFD Drag','Location','NorthWest')
xlabel(' Coefficients of Lift ')
ylabel(' Coefficient of Drag')
title('Drag Polar For the Glider')
grid minor

% - - - 

figure(5)

plot(CL_3DWing_Estimated,CD_theWholeAirplane_Polar,'*-','LineWidth',1)
hold on
% plot(CL_CFD,CD_CFD,'*-','LineWidth',1)
% hold on
% plot(CL_3DWing_Estimated,WingDrag,'*-','LineWidth',1)
% hold on
refline(0)
hold off

%legend('The whole airplane drag polar','estimated CFD Drag','3D Wing Drag','Location','NorthWest')
legend('The whole airplane drag polar','estimated CFD Drag','Location','NorthWest')
xlabel(' Coefficients of Lift ')
ylabel(' Coefficient of Drag')
title('Drag Polar Comparison')

grid minor
