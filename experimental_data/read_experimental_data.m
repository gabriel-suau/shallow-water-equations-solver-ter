clear

% Load experimental data from the .mat files.
topo = load("topo_data_17OCT.mat");
PT1 = load("data_PT1.mat");
PT2 = load("data_PT2.mat");
PT35 = load("data_sync_PT3-5.mat");

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
% We also take the first value as the new origin.

nmin = 1544631;
nmax = nmin + 9600;
ntot = nmax-nmin+1;
for i=nmin:nmax
  PT35.time(i) = [(i-nmin)/16.];
end

nmin2 = nmin;
nmax2 = nmin2 + 6000;
ntot2 = nmax2-nmin2+1;
for i=nmin2:nmax2
  PT1.time(i) = [(i-nmin)/10.];
  PT2.time(i) = [(i-nmin)/10.];
end

water_height_1 = [PT1.time(nmin2:nmax2) PT1.h_hyd(nmin2:nmax2)];
water_height_2 = [PT2.time(nmin2:nmax2) PT2.h_hyd(nmin2:nmax2)];
water_height_3 = [PT35.time(nmin:nmax) PT35.h_hyd(nmin:nmax,1)];
water_height_4 = [PT35.time(nmin:nmax) PT35.h_hyd(nmin:nmax,2)];
water_height_5 = [PT35.time(nmin:nmax) PT35.h_hyd(nmin:nmax,3)];

% We need to remove the effects of tide in experimental data.
% To do that, we make a polynomial regression to get the general trend
% and substract it to keep only the height variation due to waves.

coef_1 = polyfit(water_height_1(1:end,1),water_height_1(1:end,2),2);
coef_2 = polyfit(water_height_2(1:end,1),water_height_2(1:end,2),2);
coef_3 = polyfit(water_height_3(1:end,1),water_height_3(1:end,2),2);
coef_4 = polyfit(water_height_4(1:end,1),water_height_4(1:end,2),2);
coef_5 = polyfit(water_height_5(1:end,1),water_height_5(1:end,2),2);

regression_1 = [coef_1(1)*water_height_1(1:end,1).^2+coef_1(2)*water_height_1(1:end,1)];
regression_2 = [coef_2(1)*water_height_2(1:end,1).^2+coef_2(2)*water_height_2(1:end,1)];
regression_3 = [coef_3(1)*water_height_3(1:end,1).^2+coef_3(2)*water_height_3(1:end,1)];
regression_4 = [coef_4(1)*water_height_4(1:end,1).^2+coef_4(2)*water_height_4(1:end,1)];
regression_5 = [coef_5(1)*water_height_5(1:end,1).^2+coef_5(2)*water_height_5(1:end,1)];

water_height_1(1:end,2) -= regression_1(1:end);
water_height_2(1:end,2) -= regression_2(1:end);
water_height_3(1:end,2) -= regression_3(1:end);
water_height_4(1:end,2) -= regression_4(1:end);
water_height_5(1:end,2) -= regression_5(1:end);

% We have now the water level relative to the pressure sensors.
% The sensors are located at a height equal to delta_b from the ground,
% hence, we add this delta_b to get the real water height.

delta_b_1 = 0.15;
delta_b_2 = 0.16;
delta_b_3 = 0.14;
delta_b_4 = 0.11;
delta_b_5 = 0.12;

water_height_1(1:end,2) += delta_b_1; %- min(min(water_height));
water_height_2(1:end,2) += delta_b_2; %- min(min(water_height));
water_height_3(1:end,2) += delta_b_3; %- min(min(water_height));
water_height_4(1:end,2) += delta_b_4; %- min(min(water_height));
water_height_5(1:end,2) += delta_b_5; %- min(min(water_height));

% Write the topography and the corrected water height into csv files.
csvwrite("topography.csv",[transpose(topo.x) transpose(topo.z)]);
csvwrite("water_height_1.csv", water_height_1);
csvwrite("water_height_2.csv", water_height_2);
csvwrite("water_height_3.csv", water_height_3);
csvwrite("water_height_4.csv", water_height_4);
csvwrite("water_height_5.csv", water_height_5);
