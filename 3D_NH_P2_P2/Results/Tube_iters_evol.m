clear all; close all; clc; fclose('all'); 

%% Format plot
formats1pec = '%f';

%% Load data
fileID1 = fopen('FP_Iterations.oisd','r');
data1 = fscanf(fileID1,formats1pec);
ts1 = data1(1:2:end-3); 
iter1 = data1(2:2:end-2);

%%
figure 
plot(ts1(1:1:end-1), iter1(1:1:end-1),'or');