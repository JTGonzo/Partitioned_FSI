clear all; close all; clc;

fontSize = 15;

%fileID = fopen('Flap_Vel.othd','r');
fileID = fopen('Flap_Disp.othd','r');
formatSpec = '%f';
data = fscanf(fileID,formatSpec);

t = data(1:3:end-2); 
dx = data(2:3:end-1);
dy = data(3:3:end);
mag = sqrt(dx.^2+dy.^2);

figure(1)
grid on;
plot(t, dx,'-b','LineWidth',1);
xlabel('$tU_0/D$','Fontsize',fontSize,'Interpreter','Latex');
ylabel('$A_{x}/D$','Fontsize',fontSize,'Interpreter','Latex');

figure(2)
grid on;
plot(t, dy,'-r','LineWidth',1);
xlabel('$tU_0/D$','Fontsize',fontSize,'Interpreter','Latex');
ylabel('$A_{y}/D$','Fontsize',fontSize,'Interpreter','Latex');

figure(3)
grid on;
plot(t, mag,'-g','LineWidth',1);
xlabel('$tU_0/D$','Fontsize',fontSize,'Interpreter','Latex');
ylabel('$||\hat{d}||/D$','Fontsize',fontSize,'Interpreter','Latex');

