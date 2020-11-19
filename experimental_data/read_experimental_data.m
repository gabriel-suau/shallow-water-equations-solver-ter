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
figure(1)
plot(topo.x,topo.z)
hold on
plot(0,-1.19, "*")    #PT1
hold on
plot(34.83,0.07, "*") #PT2
hold on
plot(52.88,0.74, "*") #PT3
hold on
plot(57.90,0.83, "*") #PT4
hold on
plot(62.08,0.96, "*") #PT5

% Plot the height of water as a function of time for each sensor.
figure(2)
plot(PT1.time, PT1.h_hyd)
hold on
plot(PT2.time, PT2.h_hyd)
hold on
plot(PT35.time, PT35.h_hyd)
