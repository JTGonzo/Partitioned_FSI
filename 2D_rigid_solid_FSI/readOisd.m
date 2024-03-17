clc
clear all

F1 = fopen('rigidV3.oisd','r');
N = 2000 ;
dt = 0.1 ;
time = [dt:dt:N*dt];
U = 1 ;
D = 1 ;
rho = 1000 ;
Fct = 0.5*rho*U^2*D*1 ;

NO1 = 1500;
NO2 = N;

for i=1:N
    fgets(F1);
    fgets(F1);
    fgets(F1);
    fgets(F1);
    Force(i,:) = fscanf(F1,'%e %e',[2 1]);
    fgets(F1);
    fgets(F1);
    Disp(i,:) = fscanf(F1,'%e %e',[2 1]);
    fgets(F1);
end

figure(1)
plot(time(NO1:NO2),Force(NO1:NO2,1)./Fct,'-r',time(NO1:NO2),Force(NO1:NO2,2)./Fct,'-k');
hold on;
legend('Drag','Lift');
xlim([0 N*dt]);

figure(2)
plot(time(NO1:NO2),Disp(NO1:NO2,1),'-r',time(NO1:NO2),Disp(NO1:NO2,2),'-k');
xlim([0 N*dt]);

 
Cd = mean(Force(NO1:NO2,1))./Fct
Cl = mean(Force(NO1:NO2,2))./Fct
Clmax = max(Force(NO1:NO2,2))./Fct
Clrms = rms(Force(NO1:NO2,2) - mean(Force(NO1:NO2,2)))./Fct
Cdrms = rms(Force(NO1:NO2,1) - mean(Force(NO1:NO2,1)))./Fct
Aymax = max(Disp(NO1:NO2,2))
Axrms = rms(Disp(NO1:NO2,1)-mean(Disp(NO1:NO2,1)))

% Reference Data
% "On the vortex-induced oscillations of a freely vibrating cylinder in the
% vicinity of a stationary plane wall"
% Li Zhong, Weigang Yao, Kun Yang, Rajeev Jaiman, Boo Cheong Khoo
% Journal of Fluids and Structures 65 Pg. 495-526 (2016)
%
% Re = 200, Ur = 5, c = 0, m* = 10
%
% rho = 1000, mu = 5, m = 7853.9816, kx = ky = 12402.51062, c = 0
%
% mean Cd = 2.0551 - 2.1311
% rms Cl = 0.0893 - 0.0924
% rms Ax/D = 0.0087 - 0.0092
% max Ay/D = 0.5498 - 0.5693
