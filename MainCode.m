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
Iapp = 10^(-9);               % mA
u(603) = 1;                  % APPLYING CURRENT ARBITRARILY FOR NOW

% Transient Iapp
% Damped sinusoid
Iapp1 = @(t) (10^(-9))*sin(0.001*t)./exp(t./10000);

% Step function
Iapp2 = @(t) (10^(-9))*(t>=10000 & t<=20000);


% Steady-state voltage (proof of concept) ---------------------------------
vSS = -inv(A)*B*(u.*Iapp);

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
tspan = [0 5*10^4];

tic
[t,v] = ode23(@(t,v) A*v + B*(u.*Iapp),tspan,v0);
toc

tau = Rm.*Cm;

figure(3)
clf
plot(t./tau, v(:,1))
ylabel('Membrane potential at the soma [mV]'); xlabel('Dimensionless Time');

save('1regular.mat','t','v');


%% Transient current
% Damped sinusoid
v0 = zeros(numel(num),1);
tspan = [0 5*10^4];

tic
[t1,v1] = ode23(@(t,v) A*v + B*(u.*Iapp1(t)),tspan,v0);
toc

tau = Rm.*Cm;

figure(4)
clf
plot(t1./tau, v1(:,1))
ylabel('Membrane potential at the soma [mV] in response to damped sinusoid'); xlabel('Dimensionless Time');

save('2damped.mat','t1','v1');

%% Step current

v0 = zeros(numel(num),1);
tspan = [0 5*10^4];

tic
[t2,v2] = ode23(@(t,v) A*v + B*(u.*Iapp2(t)),tspan,v0);
toc

tau = Rm.*Cm;

figure(5)
clf
plot(t2./tau, v2(:,1))
ylabel('Membrane potential at the soma [mV] in response to step current'); xlabel('Dimensionless Time');

save('3step.mat','t2','v2');


%% Shoval calcium model, LTD

v0 = zeros(numel(num),1);
tspan = [0 5*10^4];

Iapp = [];

tic
[t3,v3] = ode23(@(t,v) A*v + B*(u.*Iapp),tspan,v0);
toc

tau = Rm.*Cm;

figure(6)
clf
plot(t3./tau, v3(:,1))
ylabel('Membrane potential at the soma [mV], Shouval model, LTD'); xlabel('Dimensionless Time');

save('4shouvalLTD.mat','t3','v3');


%% Shouval calcium model, LTP

v0 = zeros(numel(num),1);
tspan = [0 5*10^4];

Iapp = [];

tic
[t4,v4] = ode23(@(t,v) A*v + B*(u.*Iapp),tspan,v0);
toc

tau = Rm.*Cm;

figure(7)
clf
plot(t4./tau, v4(:,1))
ylabel('Membrane potential at the soma [mV], Shouval model, LTP'); xlabel('Dimensionless Time');

save('5shouvalLTP.mat','t4','v4');


%% Shouval calcium model (LTP) with NE, 0 min post-NE

v0 = zeros(numel(num),1);
tspan = [0 5*10^4];

Iapp = [];

tic
[t5,v5] = ode23(@(t,v) A*v + B*(u.*Iapp),tspan,v0);
toc

tau = Rm.*Cm;

figure(8)
clf
plot(t5./tau, v5(:,1))
ylabel('Membrane potential at the soma [mV], 0 min post-NE'); xlabel('Dimensionless Time');

save('5NE_0.mat','t5','v5');


%% Shouval calcium model (LTP) with NE, 30 min post-NE

v0 = zeros(numel(num),1);
tspan = [0 5*10^4];

Iapp = [];

tic
[t6,v6] = ode23(@(t,v) A*v + B*(u.*Iapp),tspan,v0);
toc

tau = Rm.*Cm;

figure(9)
clf
plot(t6./tau, v6(:,1))
ylabel('Membrane potential at the soma [mV], 30 min post-NE'); xlabel('Dimensionless Time');

save('5NE_0.mat','t6','v6');

%% Shouval calcium model (LTP) with NE, 90 min post-NE

v0 = zeros(numel(num),1);
tspan = [0 5*10^4];

Iapp = [];

tic
[t7,v7] = ode23(@(t,v) A*v + B*(u.*Iapp),tspan,v0);
toc

tau = Rm.*Cm;

figure(10)
clf
plot(t7./tau, v7(:,1))
ylabel('Membrane potential at the soma [mV], 90 min post-NE'); xlabel('Dimensionless Time');

save('5NE_0.mat','t7','v7');









%% UNUSED
% %% Time-varying current
% 
% % dvdt = M*D*r(t) ?
% 
% % (2.39) dydt = inv(M)*A*M*y + inv(M)*B*u 
% %        v = M*y
% tspan = [0,5e4]; % µs
% 
% M = diag(cm.^(-1/2));
% D = M\A*M;
% e = eig(D);
% capLamb = diag(e);
% % or (2.43) capLamb = inv(D)*(inv(M)*A*M)*D ?
% F = D'*inv(M)*B;
% 
% init = [];      % Initial conditions
% 
% %%
% tic
% [t1,v1] = ode23(@(t,v) M*D*r,tspan,init);
% toc



