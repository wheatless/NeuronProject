%% Compartmental Model
clear all
close all
clc

% Load dataset ------------------------------------------------------------
% num is label for that point
% type is what the point is classified as
% x, y, z are coordinates
% r is radius of point
% par is parent for that point

% All dendrites
load('dataset.mat');
data = dataset;

num  = data(:,1);
type = data(:,2);
x    = data(:,3);   % cm
y    = data(:,4);   % cm
z    = data(:,5);   % cm
r    = data(:,6);   % cm
par  = data(:,7);   % parent index

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

Ri   = 100;            % Ohm-cm
Rm   = 10000;          % Ohm-cm^2
Cm   = 1;              % muF/cm^2
Iapp = 1*10^(-9);              % mA    


% 3D Visualization --------------------------------------------------------

% All dendrites
t=0;
figure(1);
for i = num'
    hold on
    ind = find(par == i)';
    % Origin
    for j = ind
        % Color of line, based on type
        if type(j) == 1 %soma
            LineSpec = 'r-';
        elseif type(j) == 3 % basal dendrite
            LineSpec = 'g-';
        elseif type(j) == 4 % apical dendrite
            LineSpec = 'm-';
        end
        plot3([x(i), x(j)], [y(i), y(j)], [z(i), z(j)], LineSpec)
        
    end
    
end
plot3(x(1),y(1),z(1),'r.','MarkerSize',15);

% Compartment lengths and length constraint--------------------------------

% Taking true length of soma from -1 to 1 to be length between num = 1 and num = 2

% Compartment lengths in cm
l = zeros(size(num));
l(1) = sqrt((x(2)-x(1)).^2 + (y(2)-y(1)).^2 + (z(2)-z(1)).^2);
for i = 2:numel(num)
    xDist = x(i) - x(par(i));
    yDist = y(i) - y(par(i));
    zDist = z(i) - z(par(i));
    
    l(i) = sqrt(xDist^2 + yDist^2 + zDist^2);
end

% Lambdas
lamb = sqrt(((r).*Rm)./(2*Ri)); % in cm

% Electrotonic lengths
L = l./lamb;


% Defining conductances and capacitances ----------------------------------
cm = 2*pi*r.*l*Cm;
gi = (pi*r.^2)./(l*Ri);
gm = 2*pi*(r.*l)./Rm;

gi(1) = 0;


% A matrix ----------------------------------------------------------------

A = zeros(numel(num));

% Conductances in A
for i = num'
    % On (i,i)th spot
    A(i,i) = A(i,i) - gi(i) - gm(i);
    
    if i == 1
        continue
    end
    
    % On spots related to parents
    parents = par(i);
    
    for j = parents
        A(j,j) = A(j,j) - gi(i);
        A(j,i) = A(j,i) + gi(i);
        A(i,j) = A(i,j) + gi(i);
    end
    
end

% Divide by capacitance
for i = num'
    A(i,:) = A(i,:) ./ cm(i);
end


% B matrix ----------------------------------------------------------------
B = diag(1./cm);


% u matrix ----------------------------------------------------------------
u = zeros(numel(num),1);
u(603) = Iapp;                  % APPLYING CURRENT ARBITRARILY FOR NOW

% Steady-state voltage (proof of concept) ---------------------------------
vSS = -inv(A)*B*u;

% Plotting steady state voltage as func of dimensionless dist from soma ---
total = zeros(size(num));
figure(2); hold on;
for j = 2:numel(num)
    
    total(j) = L(j);
    
    k = par(j);
    
    while k ~= 1
        
        total(j) = total(j)+L(k);
        
        k = par(k);
    end
    
    plot([total(par(j)), total(j)], [vSS(par(j)) vSS(j)],'b')
    
end
hold off
xlabel('Dimensionless distance from soma'); ylabel('Voltage (mV)');
title('Steady state voltage');


%% Voltage over time (ODEs)

% Constant current
v0 = zeros(numel(num),1);


tic
[t,v] = ode23(@(t,v) A*v + B*u,tspan,v0);
toc

tau = Rm.*Cm;

figure(3)
clf
plot(t./tau, v(:,1))
ylabel('Membrane potential at the soma [mV]'); xlabel('Dimensionless Time');


%% Time-varying current

% dvdt = M*D*r(t) ?

% (2.39) dydt = inv(M)*A*M*y + inv(M)*B*u 
%        v = M*y
tspan = [0,5e4]; % µs

M = diag(cm.^(-1/2));
D = M\A*M;
e = eig(D);
capLamb = diag(e);
% or (2.43) capLamb = inv(D)*(inv(M)*A*M)*D ?
F = D'*inv(M)*B;

init = [];      % Initial conditions

%%
tic
[t1,v1] = ode23(@(t,v) M*D*r,tspan,init);
toc



