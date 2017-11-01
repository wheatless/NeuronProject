% 
clear all
close all
clc

% Constants
% Defined according to Shouval 2001
f = 100;                 % [Hz] Stimulation frequency
%tauCa = 0.001;          % Decay time constant of Ca2+ in spine
%tauf = 0.001;           % Fast: Time constant of Ca2+ current thru NMDAR
%taus = 0.001;           % Slow: Time constant of Ca2+ current thru NMDAR
Mg = 1;                  % [mM] Magnesium concentration 
Vr = 130;                % [mV] Reversal potential of Ca2+
a = 1;                   % [Hz/mV] Proportionality constant relating stim frequency to voltage
b = -100;                % [mV] Y-intercept of V(f) curve

V = @(f) a.*f + b;              % [mV] Postsynaptic potential
%Nf = func of fast NMDA rec;             % Fast: Magnitude of Ca2+ current thru NMDAR
%Ns = func of slow NMDA rec;             % Slow: Magnitude of Ca2+ current thru NMDAR



B = @(V) 1./(1 + exp(-0.062*V).*(Mg/3.57));
H = @(V) -B(V).*(V-Vr);
% Gnmda = tauCa.*(tauf.*Nf + taus.*Ns);   % Gain of Ca2+ influx thru NMDA rec. They use 
                                        % Gnmda = 0.01 and 0.03

Gnmda = 0.01;

v = V(f);                               % Voltage [mV]
CaSS = @(f) H(V(f)).*f.*Gnmda;          % Ca2+ flux [mM/s]?

% H(v) is in Volt/mM
% f is in 1/s
% Gnmda is in s^2 ? What is Nf and Ns in? Gnmda should be unitless?
% Thus [Ca2+] is in Volt*s/mM ??

% Plots

% Voltage vs. Freq
freq = 0:100;
volt = V(freq);
figure;
plot(freq,volt); xlabel('Freqency [Hz]'); ylabel('Voltage [mV]');
title('Voltage as function of stimulation frequency');

% H as function of V
v1 = -90:75;
figure;
plot(v1,H(v1)); xlabel('Voltage [mV]'); ylabel('H(v)');
title('H as a function of Voltage');

% B as a function of V
figure;
v2 = -100:75;
plot(v2,B(v2)); xlabel('Voltage [mV]'); ylabel('B');
title('B as a function of Voltage');

% Ca2+ Steady-State as function of stimulation frequency
f1 = 0:300;
figure;
plot(f1,CaSS(f1)); xlabel('Frequency [Hz]'); ylabel('[Ca^2^+] [mM]');
title('Steady state [Ca^2^+] as function of stimulation frequency');

%% Next step

%Functions for phosphotase and kinase rates
EK = @(ca) 1+100.*(ca.^2)./(64 + (ca.^2));
EP = @(ca) 1+30.*(ca.^2)./(1+(ca.^2));

%Plot rates
ca = 0:30;
EP1 = EP(ca);
EP2 = EP1;
EK1 = EK(ca);
EK2 = EK1;
figure(1);
plot(ca, EP1, ca, EK1);
legend('K1&K2', 'P1&P2', 'Location', 'Best');
xlabel('Calcium'); ylabel('Rate (1/s)');
title('Enzymatic Activity (no NE)');

%Calculate and plot fraction in each state
A = (EP1.*EP2)./((EK2+EP2).*(EK1+EP1));
Ap1 = EK1.*EP2./((EK2+EP2).*(EK1+EP1));
Ap2 = EK2.*EP1./((EK2+EP2).*(EK1+EP1));
Ap1p2 = (EK1.*EK2)./((EK2+EP2).*(EK1+EP1));

figure(2); hold on;
plot(ca, A, ca, Ap1, ca,Ap2, ca, Ap1p2);
legend('A', 'Ap1', 'Ap2', 'Ap1p2', 'Location', 'Best');
xlabel('Calcium'); ylabel('GLuR1 State Fractions');
title('Phosophorylation (no NE)');

%Model ampa conducatance (matches experimental)
gampa = A + 2*(Ap1+Ap2) + 4*(Ap1p2);

figure(3); 
plot(ca, gampa);
xlabel('Calcium'); ylabel('Conductance');
title('AMPA Conductance');

%% Increase norepinephrine

EK2 = EK2.*2; %Activation of CaMKII increases rate of phosphorylation of ser831
EP1 = EP1./2; %Activation of PKA increases mean open time - decrease EP1

%Recalculate fraction in states
A = (EP1.*EP2)./((EK2+EP2).*(EK1+EP1));
Ap1 = EK1.*EP2./((EK2+EP2).*(EK1+EP1));
Ap2 = EK2.*EP1./((EK2+EP2).*(EK1+EP1));
Ap1p2 = (EK1.*EK2)./((EK2+EP2).*(EK1+EP1));

figure(4); hold on;
plot(ca, A, ca, Ap1, ca,Ap2, ca, Ap1p2);
legend('A', 'Ap1', 'Ap2', 'Ap1p2', 'Location', 'Best');
xlabel('Calcium'); ylabel('GLuR1 State Fractions');
title('Phosophorylation (NE applied)');

%recalculate ampa conducatance
gampa = A + 2*(Ap1+Ap2) + 4*(Ap1p2);

%plot augmented ampa conductance 
figure(3); hold on;
plot(ca, gampa);
legend('No NE', 'NE applied', 'Location', 'Best');
