%% Compartmental Model
clear all
close all
clc

% Load dataset

% All dendrites
load('dataset.mat');
dataAll = dataset;

numAll  = dataAll(:,1);
typeAll = dataAll(:,2);
xAll    = dataAll(:,3);   % cm
yAll    = dataAll(:,4);   % cm
zAll    = dataAll(:,5);   % cm
rAll    = dataAll(:,6);   % cm
parAll  = dataAll(:,7);   % parent index

% Apical dendrites
load('datasetApical.mat');
dataApi = datasetApical;

numApi  = dataApi(:,1);
typeApi = dataApi(:,2);
xApi    = dataApi(:,3);   % cm
yApi    = dataApi(:,4);   % cm
zApi    = dataApi(:,5);   % cm
rApi    = dataApi(:,6);   % cm
parApi  = dataApi(:,7);   % parent index

%Basal dendrites
load('datasetBasal.mat');
dataBas = datasetBasal;

numBas  = dataBas(:,1);
typeBas = dataBas(:,2);
xBas    = dataBas(:,3);   % cm
yBas    = dataBas(:,4);   % cm
zBas    = dataBas(:,5);   % cm
rBas    = dataBas(:,6);   % cm
parBas  = dataBas(:,7);   % parent index

% Constants

Ri   = 1;              % Ohm-cm
Rm   = 1;              % Ohm-cm^2
Cm   = 1;              % muF/cm^2
Iapp = 1;              % mA    


% 3D Visualization

% All dendrites
for i = numAll
    ind = find(parAll==i);
    
    plot3([xAll(i), xAll(ind)], [yAll(i), yAll(ind)], [zAll(i), zAll(ind)]);
    
end
hold on
% Soma
plot3(xAll(1),yAll(1),zAll(1),'r.','MarkerSize',10);








