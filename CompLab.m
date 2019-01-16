%% info


%% Housekeeping

clear;
clc;
close all;

%% Read data

%CFD Data
TempsetCFD = xlsread("Tempest UAS CFD flight data CFD.xlsx");
alphaCFD = TempsetCFD(:,1);
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

% Alpha when slope is == 0 
Alpha0 = mean([Alpha2D(Postiv);Alpha2D(Postiv-1)]);

% a0 = cl / alpha2D

% pick point in the middle

[ r c ] = size(Cl_2D);
MidPoint = floor(r/2);


a0 = (Cl_2D(MidPoint+1)-Cl_2D(MidPoint))/((Alpha2D(MidPoint+1)-Alpha2D(MidPoint)));

% the lift curve slope for 3D

a3D_liftCurveSlope = (a0)/(1+ ( ( 57.3 * a0 ) / ( pi * EfficRatio * AspectRatio) )) ;

CL_3D_Estimated = a3D_liftCurveSlope .* ( Alpha2D - Alpha0);

%% Induced drag

InducedDrag = (CL_3D_Estimated).^2 ./ (pi*EfficRatio*AspectRatio);

% the wing drag 
WingDrag = Cd_2D + InducedDrag;
AlphaMin = Alpha2D(find(WingDrag == min(WingDrag)));

% Oswald efficiency e0

e0 = 1.78 * ( 1 - 0.045 .* AspectRatio.^(0.68)) - 0.64 ;
k1 = 1 / ( pi * e0 * AspectRatio ) ;

%% plot

figure(1)

plot(Alpha2D,CL_3D_Estimated)
hold on
plot(alphaCFD,CL_CFD)
hold on
plot(Alpha2D,Cl_2D)
hold off

legend('3D Calculated','CFD','2D')
