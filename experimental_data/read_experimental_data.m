clear

% Load experimental data from the .mat files.
topo = load("topo_data_17OCT.mat");
PT1 = load("data_PT1.mat");
PT2 = load("data_PT2.mat");
PT35 = load("data_sync_PT3-5.mat");

% Get time vector.
PT1.time = PT1.time - 738081;
PT2.time = PT2.time - 738081;
PT35.time = PT35.time - 738081;

% Plot the topography and the position of the sensors.

%figure(1)
%plot(topo.x,topo.z)
%hold on
%plot(0,-1.19, "*")    #PT1
%hold on
%plot(34.83,0.07, "*") #PT2
%hold on
%plot(52.88,0.74, "*") #PT3
%hold on
%plot(57.90,0.83, "*") #PT4
%hold on
%plot(62.08,0.96, "*") #PT5
%
% Plot the height of water as a function of time for each sensor.
%figure(2)
%plot(PT1.time, PT1.h_hyd)
%hold on
%plot(PT2.time, PT2.h_hyd)
%hold on
%plot(PT35.time, PT35.h_hyd)

% In our code, we are going to work with only a small portion of the data.
% We chose the time period [0.56027, 0.62031] and get the corresponding
% indices : 0.56027 -> PT35(113000) ; 0.62031 -> PT35(196000)

nmin = 113000;
nmax = 196000;
ntot=nmax-nmin+1;
water_height = [PT35.time(nmin:nmax) PT35.h_hyd(nmin:nmax,1) PT35.h_hyd(nmin:nmax,2) PT35.h_hyd(nmin:nmax,3)];

% We need to remove the effects of tide in experimental data.
% To do that, we make a polynomial regression to get the general trend
% and substract it to keep only the height variation due to waves.

coef_3 = polyfit(water_height(1:end,1),water_height(1:end,2),2);
coef_4 = polyfit(water_height(1:end,1),water_height(1:end,3),2);
coef_5 = polyfit(water_height(1:end,1),water_height(1:end,4),2);

regression_3 = [coef_3(1)*water_height(1:end,1).^2+coef_3(2)*water_height(1:end,1)+coef_3(3)];
regression_4 = [coef_4(1)*water_height(1:end,1).^2+coef_4(2)*water_height(1:end,1)+coef_4(3)];
regression_5 = [coef_5(1)*water_height(1:end,1).^2+coef_5(2)*water_height(1:end,1)+coef_5(3)];

water_height(1:end,2) -= regression_3(1:end);
water_height(1:end,3) -= regression_4(1:end);
water_height(1:end,4) -= regression_5(1:end);

% We have now the water level relative to the pressure sensors.
% The sensors are located at a height equal to delta_b from the ground,
% hence, we add this delta_b to get the real water height.

delta_b_3 = 0.14;
delta_b_4 = 0.11;
delta_b_5 = 0.12;

water_height(1:end,2) += delta_b_3 - min(min(water_height));
water_height(1:end,3) += delta_b_4 - min(min(water_height));
water_height(1:end,4) += delta_b_5 - min(min(water_height));

% Write the topography and the corrected water height into csv files.
csvwrite("topography.csv",[transpose(topo.x) transpose(topo.z)]);
csvwrite("water_height_35.csv", water_height);
