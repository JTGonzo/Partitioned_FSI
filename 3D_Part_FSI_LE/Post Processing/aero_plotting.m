clear all; close all; clc;

fontSize = 15;

% fileID = fopen('AeroForces.txt','r');
fileID = fopen('AerodynamicForcesFSI_Final.txt','r');
formatSpec = '%f';
data = fscanf(fileID,formatSpec);

t = data(1:4:end-3); 
Fx = data(2:4:end-2);
Fy = data(3:4:end-1);
Fz = data(4:4:end);
%mag = sqrt(Fx.^2+Fy.^2+Fz.^2);

figure(1)
grid on;
plot(t, Fx,'-b','LineWidth',1);
xlabel('$t$','Fontsize',fontSize,'Interpreter','Latex');
ylabel('$F_{x}$','Fontsize',fontSize,'Interpreter','Latex');

figure(2)
grid on;
plot(t, Fy,'-r','LineWidth',1);
xlabel('$t$','Fontsize',fontSize,'Interpreter','Latex');
ylabel('$F_{y}$','Fontsize',fontSize,'Interpreter','Latex');

figure(3)
grid on;
plot(t, Fz,'-m','LineWidth',1);
xlabel('$t$','Fontsize',fontSize,'Interpreter','Latex');
ylabel('$F_{z}$','Fontsize',fontSize,'Interpreter','Latex');

% figure(4)
% grid on;
% plot(t, mag,'-g','LineWidth',1);
% xlabel('$t$','Fontsize',fontSize,'Interpreter','Latex');
% ylabel('$||\hat{d}||$','Fontsize',fontSize,'Interpreter','Latex');

