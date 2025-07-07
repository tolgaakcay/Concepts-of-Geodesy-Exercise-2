clc
clear 
format long
% Loading the data into MATLAB
load coast
load data
load lageos
load etalon
lageos = table2array(lageos);
lageos = str2double(lageos);
etalon = table2array(etalon);
etalon = str2double(etalon);
% Lageos Latitude/Longitude Calculations
xL = lageos(1:720,8);
yL = lageos(1:720,9);
zL = lageos(1:720,10);
longL = atan2(yL,xL);
longL = 180*longL/pi;
latL = atand(zL./(sqrt(yL.^2 + xL.^2)));
% Etalon Latitude/Longitude Calculations
xE = etalon(1:96,8);
yE = etalon(1:96,9);
zE = etalon(1:96,10);
longE = atan2(yE,xE);
longE = 180*longE/pi;
latE = atand(zE./sqrt(yE.^2 + xE.^2));
% Plotting the first day's ground track for the satellites
figure(1)
plot(longL,latL,'.');
hold on
plot(long,lat);
legend('Location of Lageos1');
xlabel('Longitude');
ylabel('Latitude');
title('One Day Track Map for Lageos 1 (05.09.2021)')
grid on

figure(2)
plot(longE,latE,'.');
hold on
plot(long,lat);
legend('Location of Etalon1');
xlabel('Longitude');
ylabel('Latitude');
title('One Day Track Map for Etalon 1 (05.09.2021)')
grid on 

% Finding UTM Zone for each satellite's first position
zoneL = utmzone(latL(1),longL(1));
zoneE = utmzone(latE(1),longE(1));

% Getting Polar coordinates Xp and Yp from Ex1's dataset
dset = table2array(data);
pos = find((dset(:,1) == 2021) & (dset(:, 2) == 09) & (dset(:, 3) == 05));
Xp = dset(pos, 5);
Yp = dset(pos, 6);
% Calculating the MJD(UTC), UTC time and UT1 time
MJD_1 = dset(pos, 4);
MJD_UTC_L = [MJD_1]; MJD_UTC_E = [MJD_1];
for i = 1:719
     MJD_UTC_L(end+1) = MJD_UTC_L(i) + 0.00138888889;
end
for i = 1:95
     MJD_UTC_E(end+1) = MJD_UTC_E(i) + 0.0104166667;
end
UTC_L = {}; UTC_E = {};
UT1_L = {}; UT1_E = {};
DUT = dset(pos,7);
for i = 1:720
    UTC_L{end+1} = {lageos(i,1), lageos(i,2), lageos(i,3), lageos(i,4), lageos(i,5), lageos(i,6)};
    UT1_L{end+1} = {lageos(i,1), lageos(i,2), lageos(i,3), lageos(i,4), lageos(i,5), (lageos(i,6) + DUT)};
end
for i = 1:96
    UTC_E{end+1} = {etalon(i,1), etalon(i,2), etalon(i,3), etalon(i,4), etalon(i,5), etalon(i,6)};
    UT1_E{end+1} = {etalon(i,1), etalon(i,2), etalon(i,3), etalon(i,4), etalon(i,5), (etalon(i,6) + DUT)};
end
MJD_TT_L = MJD_UTC_L + 0.000800740741; % 69.184 seconds in days
MJD_TT_E = MJD_UTC_E + 0.000800740741;
% Finding Tu and t, then deriving s prime and ERA for the satellites
% For Lageos
JD_UTC_L = MJD_UTC_L + 2400000.5;
JD_UT1_L = JD_UTC_L + DUT;
JD_TT_L = MJD_TT_L + 2400000.5;
TuL = JD_UT1_L - 2451545;
tL = (JD_TT_L - 2451545)/36525;

sP_L = -47*tL;
ERA_L = 2*pi*(0.7790572732640 + (1.00273781191135448*TuL));

% For Etalon
JD_UTC_E = MJD_UTC_E + 2400000.5;
JD_UT1_E = JD_UTC_E + DUT;
JD_TT_E = MJD_TT_E + 2400000.5;
TuE = JD_UT1_E - 2451545;
tE = (JD_TT_E - 2451545)/36525;

sP_E = -47*tE;
ERA_E = 2*pi*(0.7790572732640 + (1.00273781191135448*TuE));

% Rotation Matrixes
a = 0;
R1 = [1,0,0;0,cosd(a),sind(a);0,(-sind(a)),cosd(a)];
R2 = [cosd(a), 0, (-sind(a)); 0, 1, 0; sind(a), 0, cosd(a)];
R3 = [cosd(a), sind(a), 0; (-sind(a)), cosd(a), 0; 0, 0, 1];
% R - Earth Rotation Matrix
R_L = {cosd(-(ERA_L)), sind(-(ERA_L)), 0; (-sind(-(ERA_L))), cosd(-(ERA_L)), 0; 0, 0, 1};
R_E = {cosd(-(ERA_E)), sind(-(ERA_E)), 0; (-sind(-(ERA_E))), cosd(-(ERA_E)), 0; 0, 0, 1};

% W - Polar Motion Matrix
R1_L = [1,0,0;0,cosd(Yp),sind(Yp);0,(-sind(Yp)),cosd(Yp)];
R2_L = [cosd(Xp), 0, (-sind(Xp)); 0, 1, 0; sind(Xp), 0, cosd(Xp)];
R3_L = {cosd(-(sP_L)), sind(-(sP_L)), 0; (-sind(-(sP_L))), cosd(-(sP_L)), 0; 0, 0, 1};
W_L = R1_L.*R2_L*R3_L;


