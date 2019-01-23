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


Airfoil2D_Data = xlsread("Airfoil2D Data.xlsx");
Alpha2D = Airfoil2D_Data(:,1);
Cl_2D = Airfoil2D_Data(:,2);
Cd_2D = Airfoil2D_Data(:,3);


AspectRatio = 16.5 ; 
EfficRatio = 0.9;

%% find where lift is == 0 

% for this we will take the avg between the negative and possitive

% first positive number
Postiv = find(Cl_2D>0,1);

% Alpha when lift is == 0 
Alpha0 = mean([Alpha2D(Postiv);Alpha2D(Postiv-1)]);
%Alpha0 = (Cl_2D(Postiv) - Cl_2D(Postiv - 1))/(Alpha2D(Postiv) - Alpha2D(Postiv - 1)) * (Alpha2D(Postiv));

% a0 = cl / alpha2D

% pick point in the middle

[ r c ] = size(Cl_2D);
MidPoint = floor(r/2);


a0 = (Cl_2D(MidPoint+1)-Cl_2D(MidPoint))/((Alpha2D(MidPoint+1)-Alpha2D(MidPoint)));

% the lift curve slope for 3D

a3D_liftCurveSlope = (a0)/(1+ ( ( 57.3 * a0 ) / ( pi * EfficRatio * AspectRatio) )) ;

CL_3D_Estimated = a3D_liftCurveSlope .* ( Alpha2D - Alpha0);

%% Induced drag, and  wing drag

InducedDrag = (CL_3D_Estimated).^2 ./ (pi*EfficRatio*AspectRatio);

% the wing drag 
WingDrag = Cd_2D + InducedDrag;

%% the whole aircraft drag

% find where the the minumm drage happens on the whole wing, get CL
% corresponding. 
CL_Min_D = CL_3D_Estimated(find(WingDrag == min(WingDrag)));


% Oswald efficiency e0

e0 = 1.78 * ( 1 - 0.045 .* AspectRatio.^(0.68)) - 0.64 ;
k1 = 1 / ( pi * e0 * AspectRatio ) ;


% total polar drag: for the whole aircraft.

CD_Polar = min(WingDrag) + k1.*((CL_3D_Estimated - CL_Min_D).^2) ; 


%% L/D : is it 3d just wing or the whoel aircraft/


L_D_WholeAirplane = CL_3D_Estimated ./ WingDrag ; 
L_D_CFD = (CL_CFD./CD_CFD);
%% plot

figure(1)

plot(Alpha2D,CL_3D_Estimated,'-.','LineWidth',2.5)
hold on
plot(AlphaCFD,CL_CFD,'-.','LineWidth',2.5)
hold on
plot(Alpha2D,Cl_2D,'-.','LineWidth',2.5)
hold on
refline(0)
hold off

legend('3D Calculated','CFD','2D','Location','SouthEast')
xlabel(' \alpha \circ ')
ylabel(' Coefficient of Lift')
title('Lift Curve Comparison')
grid minor


% - - -


figure(2)

plot(CL_3D_Estimated,WingDrag,'-.','LineWidth',2.5)
hold on
plot(CL_3D_Estimated,CD_Polar,'-.','LineWidth',2.5)
hold on
plot(CL_CFD,CD_CFD,'-.','LineWidth',2.5)
hold on
refline(0)
hold off

legend('Finite wing drag','Polar drag (the whole aircraft)','CFD Drag','Location','NorthWest')
xlabel(' Coefficients of Lift ')
ylabel(' Coefficient of Drag')
title('Drag Polar Comparison')
grid minor

% - - - 

figure(3)

plot(Alpha2D,L_D_WholeAirplane,'-.','LineWidth',2.5)
hold on
plot(AlphaCFD,L_D_CFD,'-.','LineWidth',2.5)
hold on
refline(0)
hold off

legend('3D Wings full L/D','CFD L/D','Location','NorthWest')
xlabel(' \alpha \circ ')
ylabel(' L/D')
title('L/D Comparison')
grid minor

% - - - 