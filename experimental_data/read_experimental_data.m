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

% 0.56027 -> PT35(113000) - 0.62031 ->PT35(196000)

nmin = 113000;
nmax = 196000;
valeurs = [PT35.time(nmin:nmax) PT35.h_hyd(nmin:nmax,1) PT35.h_hyd(nmin:nmax,2) PT35.h_hyd(nmin:nmax,3)];

ntot=nmax-nmin+1;
moy = mean(valeurs(1:end,2))
dt = (valeurs(ntot,1)-valeurs(1,1))
droite = [];

coef = polyfit(valeurs(1:end,1),valeurs(1:end,2),2);

droite = [coef(1)*valeurs(1:end,1).^2+coef(2)*valeurs(1:end,1)+coef(3)];
new_valeur = valeurs(1:end,2)-droite(1:end);
new_valeur2 = new_valeur-min(new_valeur)


plot(new_valeur2)
















csvwrite('valeurs.csv', new_valeur2);
