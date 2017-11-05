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
Iapp = 1.3*10^(-9);              % mA    


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
plot3(x(570),y(570),z(570),'bx','MarkerSize',10);

% Find where branches end -------------------------------------------------

ends = ~ismember(num,par);
endIndex = find(ends==1);

% Choosing end to apply current




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
Iapp = 1.3*10^(-9);               % mA
u(570) = 1;                  % APPLYING CURRENT ARBITRARILY FOR NOW

% Transient Iapp
% Damped sinusoid
Iapp1 = @(t) (1.3*10^(-9))*sin(0.001*t)./exp(t./10000);

% Step function
Iapp2 = @(t) (1.3*10^(-9))*(t>=10000 & t<=20000);


% Steady-state voltage ----------------------------------------------------
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


%% VOLTAGE OVER TIME (ODEs)

% Constant current
v0 = zeros(numel(num),1);
tspan = [0 5*10^4];

Iapp = 1.3*10^(-9);         % [mA]

tic
[t,v] = ode23(@(t,v) A*v + B*(u.*Iapp),tspan,v0);
toc

tau = Rm.*Cm;

figure(3)
clf
plot(t./tau, v(:,1))
ylabel('Membrane potential at the soma [mV]'); xlabel('Dimensionless Time');
title('Membrane potential at the soma [mV]');

save('1regulart.mat','t');
save('1regularvSS.mat','vSS');
v01 = v(:,1:400);
save('1regularv1.mat','v01');
v02 = v(:,401:800);
save('1regularv2.mat','v02');
v03 = v(:,801:1200);
save('1regularv3.mat','v03');
v04 = v(:,1201:end);
save('1regularv4.mat','v04');


load handel
sound(y,Fs)

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
title('Membrane potential at the soma [mV] in response to damped sinusoid'); xlabel('Dimensionless Time');
ylabel('Membrane potential at the soma [mV]');

save('2dampedt.mat','t1');
clear v01 v02 v03 v04
v01 = v1(:,1:400);
save('2dampedv1.mat','v01');
v02 = v1(:,401:800);
save('2dampedv2.mat','v02');
v03 = v1(:,801:1200);
save('2dampedv3.mat','v03');
v04 = v1(:,1201:end);
save('2dampedv4.mat','v04');

load handel
sound(y,Fs)

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
title('Membrane potential at the soma [mV] in response to step current'); xlabel('Dimensionless Time');
ylabel('Membrane potential at the soma [mV]');

save('3stept.mat','t2');
clear v01 v02 v03 v04
v01 = v2(:,1:400);
save('3stepv1.mat','v01');
v02 = v2(:,401:800);
save('3stepv2.mat','v02');
v03 = v2(:,801:1200);
save('3stepv3.mat','v03');
v04 = v2(:,1201:end);
save('3stepv4.mat','v04');

load handel
sound(y,Fs)

%% Channel dynamics at synapse, constant current again

% Shoval calcium model, LTD

v0 = zeros(numel(num),1);
tspan = [0 5*10^4];


Iapp = 1.3*10^(-9)*0.5887;

tic
[t3,v3] = ode23(@(t,v) A*v + B*(u.*Iapp),tspan,v0);
toc

% Steady-state voltage
clear('vSS');
vSS = -inv(A)*B*(u.*Iapp);

tau = Rm.*Cm;

figure(6)
clf
plot(t3./tau, v3(:,1))
title('Membrane potential at the soma [mV], Shouval model, LTD'); xlabel('Dimensionless Time');
ylabel('Membrane potential at the soma [mV]');

save('4shouvalLTDt.mat','t3');
save('4shouvalLTDvSS.mat','vSS');
clear v01 v02 v03 v04
v01 = v3(:,1:400);
save('4shouvalLTDv1.mat','v01');
v02 = v3(:,401:800);
save('4shouvalLTDv2.mat','v02');
v03 = v3(:,801:1200);
save('4shouvalLTDv3.mat','v03');
v04 = v3(:,1201:end);
save('4shouvalLTDv4.mat','v04');

% Plotting steady state voltage as func of dimensionless dist from soma ---
total = zeros(size(num));
figure(11); hold on;
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
title('Steady state voltage, Shouval model, LTD');


load handel
sound(y,Fs)

%% Shouval calcium model, LTP

v0 = zeros(numel(num),1);
tspan = [0 5*10^4];

Iapp = 1.3*10^(-9)*1.2729;

tic
[t4,v4] = ode23(@(t,v) A*v + B*(u.*Iapp),tspan,v0);
toc

% Steady-state voltage
clear('vSS');
vSS = -inv(A)*B*(u.*Iapp);

tau = Rm.*Cm;

figure(7)
clf
plot(t4./tau, v4(:,1))
title('Membrane potential at the soma [mV], Shouval model, LTP'); xlabel('Dimensionless Time');
ylabel('Membrane potential at the soma [mV]');

save('5shouvalLTPt.mat','t4');
save('5shouvalLTPvSS.mat','vSS')
clear v01 v02 v03 v04
v01 = v4(:,1:400);
save('5shouvalLTPv1.mat','v01');
v02 = v4(:,401:800);
save('5shouvalLTPv2.mat','v02');
v03 = v4(:,801:1200);
save('5shouvalLTPv3.mat','v03');
v04 = v4(:,1201:end);
save('5shouvalLTPv4.mat','v04');


% Plotting steady state voltage as func of dimensionless dist from soma ---
total = zeros(size(num));
figure(12); hold on;
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
title('Steady state voltage, Shouval model, LTP');

% SAVE V4 TOO
load handel
sound(y,Fs)

%% Shouval calcium model (LTP) with NE, 0 min post-NE

v0 = zeros(numel(num),1);
tspan = [0 5*10^4];

Iapp = @(t) 1.803*(1.3*10^(-9))*(t>=10000 & t<=20000);

tic
[t5,v5] = ode23(@(t,v) A*v + B*(u.*Iapp(t)),tspan,v0);
toc

tau = Rm.*Cm;

figure(8)
clf
plot(t5./tau, v5(:,1))
title('Membrane potential at the soma [mV], 0 min post-NE'); xlabel('Dimensionless Time');
ylabel('Membrane potential at the soma [mV]');


save('6NE_0.mat','t5');
clear v01 v02 v03 v04
v01 = v5(:,1:400);
save('6NE_0v1.mat','v01');
v02 = v5(:,401:800);
save('6NE_0v2.mat','v02');
v03 = v5(:,801:1200);
save('6NE_0v3.mat','v03');
v04 = v5(:,1201:end);
save('6NE_0v4.mat','v04');

load handel
sound(y,Fs)

%% Shouval calcium model (LTP) with NE, 30 min post-NE

v0 = zeros(numel(num),1);
tspan = [0 5*10^4];

Iapp = @(t) 1.6941*(1.3*10^(-9))*(t>=10000 & t<=20000);

tic
[t6,v6] = ode23(@(t,v) A*v + B*(u.*Iapp(t)),tspan,v0);
toc

tau = Rm.*Cm;

figure(9)
clf
plot(t6./tau, v6(:,1))
title('Membrane potential at the soma [mV], 30 min post-NE'); xlabel('Dimensionless Time');
ylabel('Membrane potential at the soma [mV]');

save('7NE_30.mat','t6');
clear v01 v02 v03 v04
v01 = v6(:,1:400);
save('7NE_30v1.mat','v01');
v02 = v6(:,401:800);
save('7NE_30v2.mat','v02');
v03 = v6(:,801:1200);
save('7NE_30v3.mat','v03');
v04 = v6(:,1201:end);
save('7NE_30v4.mat','v04');

load handel
sound(y,Fs)

%% Shouval calcium model (LTP) with NE, 90 min post-NE

v0 = zeros(numel(num),1);
tspan = [0 5*10^4];

Iapp = @(t) 1.3611*(1.3*10^(-9))*(t>=10000 & t<=15000);

tic
[t7,v7] = ode23(@(t,v) A*v + B*(u.*Iapp(t)),tspan,v0);
toc

tau = Rm.*Cm;

figure(10)
clf
plot(t7./tau, v7(:,1))
title('Membrane potential at the soma [mV], 90 min post-NE'); xlabel('Dimensionless Time');
ylabel('Membrane potential at the soma [mV]');

save('8NE_90.mat','t7');
clear v01 v02 v03 v04
v01 = v7(:,1:400);
save('8NE_90v1.mat','v01');
v02 = v7(:,401:800);
save('8NE_90v2.mat','v02');
v03 = v7(:,801:1200);
save('8NE_90v3.mat','v03');
v04 = v7(:,1201:end);
save('8NE_90v4.mat','v04');


load handel
sound(y,Fs)

%% Shouval calcium model (LTP), No NE

v0 = zeros(numel(num),1);
tspan = [0 5*10^4];

Iapp = @(t) (1.3*10^(-9))*(t>=10000 & t<=15000);

tic
[t8,v8] = ode23(@(t,v) A*v + B*(u.*Iapp(t)),tspan,v0);
toc

tau = Rm.*Cm;

figure;
clf
plot(t8./tau, v8(:,1))
title('Membrane potential at the soma [mV], No NE'); xlabel('Dimensionless Time');
ylabel('Membrane potential at the soma [mV]');

save('9NE_No.mat','t8');
clear v01 v02 v03 v04
v01 = v8(:,1:400);
save('9NE_Nov1.mat','v01');
v02 = v8(:,401:800);
save('9NE_Nov2.mat','v02');
v03 = v8(:,801:1200);
save('9NE_Nov3.mat','v03');
v04 = v8(:,1201:end);
save('8NE_90v4.mat','v04');


load handel
sound(y,Fs)
