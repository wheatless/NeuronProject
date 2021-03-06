% 
clear all
close all
clc

% Constants
% Defined according to Shouval 2001
f = 20;                  % [Hz] Stimulation frequency (0.5-3Hz = LTD, >10 = LTP)
fLTP = 30; 
fLTD = 2.5;

Mg = 1;                  % [mM] Magnesium concentration 
Vr = 130;                % [mV] Reversal potential of Ca2+
a = 1;                   % [Hz/mV] Proportionality constant relating stim frequency to voltage
b = -70;                 % [mV] Y-intercept of V(f) curve
Gnmda = 0.01;            % Taken from Shouval 2001


% UNUSED CONSTANTS
% tauCa = 0.001;                          % Decay time constant of Ca2+ in spine
% tauf = 0.001;                           % Fast: Time constant of Ca2+ current thru NMDAR
% taus = 0.001;                           % Slow: Time constant of Ca2+ current thru NMDAR
% Nf = func of fast NMDA rec;             % Fast: Magnitude of Ca2+ current thru NMDAR
% Ns = func of slow NMDA rec;             % Slow: Magnitude of Ca2+ current thru NMDAR
% Gnmda = tauCa.*(tauf.*Nf + taus.*Ns);   % Gain of Ca2+ influx thru NMDA rec. They use 
                                          % Gnmda = 0.01 and 0.03

% Equations
V = @(f) a.*f + b;                      % [mV] Postsynaptic potential

B = @(V) 1./(1 + exp(-0.062*V).*(Mg/3.57)); % Ca current (thru NMDAR) dependence on Mg block
H = @(V) -B(V).*(V-Vr);                     % Voltage dependence of Ca current thru NMDAR

% Values
v = V(f);                               % Voltage [mV]
CaSS = @(f) H(V(f)).*f.*Gnmda;          % Steady state [Ca2+] [mM]

Ca = CaSS(f);

CaLTP = CaSS(fLTP);
CaLTD = CaSS(fLTD);


% Plots

% Voltage vs. Freq, Shouval 2001 Fig 3a
freq = 0:100;
volt = V(freq);
figure;
plot(freq,volt); xlabel('Freqency [Hz]'); ylabel('Voltage [mV]');
title('Voltage as function of stimulation frequency');

% H as function of V, Shouval 2002 Fig 5b
v1 = -90:75;
figure;
plot(v1,H(v1)); xlabel('Voltage [mV]'); ylabel('H(v)');
title('H as a function of Voltage');
% Looks like Shouval 2002

% B as a function of V, Shouval 2002 Fig 5a
figure;
v2 = -100:75;
plot(v2,B(v2)); xlabel('Voltage [mV]'); ylabel('B');
title('B as a function of Voltage');
% Looks like Shouval 2002

% Ca2+ Steady-State as function of stimulation frequency
f1 = 0:60;
figure;
plot(f1,CaSS(f1)); xlabel('Frequency [Hz]'); ylabel('[Ca^2^+] [mM]');
title('Steady state [Ca^2^+] as function of stimulation frequency');
% Note Shouval [Ca2+] values don't go above ~40 mM

%% Phosphatases, kinases, and their effect on AMPAR conductance

%Functions for phosphotase and kinase rates
EK = @(ca) 1+100.*(ca.^2)./(64 + (ca.^2));
EP = @(ca) 1+30.*(ca.^2)./(1+(ca.^2));

%Plot rates
%ca = linspace(0,30,100);
% ca = CaSS(f1);
ca = [CaSS(0) CaSS(30)];
EP1 = EP(ca);
EP2 = EP1;
EK1 = EK(ca);
EK2 = EK1;
figure;
plot(ca, EP1, ca, EK1);
legend('P1&P2', 'K1&K2', 'Location', 'Best');
xlabel('Calcium'); ylabel('Rate (1/s)');
title('Enzymatic Activity (no NE)');
xlim([0,30]);

%Calculate and plot fraction in each state
A = (EP1.*EP2)./((EK2+EP2).*(EK1+EP1));
Ap1 = (EK1.*EP2)./((EK2+EP2).*(EK1+EP1));
Ap2 = (EK2.*EP1)./((EK2+EP2).*(EK1+EP1));
Ap1p2 = (EK1.*EK2)./((EK2+EP2).*(EK1+EP1));

figure; hold on;
plot(ca, A, ca, Ap1, ca,Ap2, ca, Ap1p2);
legend('A', 'Ap1', 'Ap2', 'Ap1p2', 'Location', 'Best');
xlabel('Calcium'); ylabel('GLuR1 State Fractions');
title('Phosophorylation (no NE)');
xlim([0,30]);

%Model ampa conducatance (matches experimental)
gampa1 = A + 2*(Ap1+Ap2) + 4*(Ap1p2);

figure; 
plot(ca, gampa1);
xlabel('Calcium'); ylabel('Conductance');
title('AMPA Conductance');
xlim([0,30]);

% AMPA conductance vs frequency

figure
plot(f1,gampa1/gampa1(1));
xlabel('Frequency [Hz]'); ylabel('Conductance');
title('AMPA Conductance vs. Frequency');

%% Increase norepinephrine

EK2 = EK2.*100; %Activation of CaMKII increases rate of phosphorylation of ser831
EP1 = EP1./100; %Activation of PKA increases mean open time - decrease EP1

%Recalculate fraction in states
A = (EP1.*EP2)./((EK2+EP2).*(EK1+EP1));
Ap1 = EK1.*EP2./((EK2+EP2).*(EK1+EP1));
Ap2 = EK2.*EP1./((EK2+EP2).*(EK1+EP1));
Ap1p2 = (EK1.*EK2)./((EK2+EP2).*(EK1+EP1));

figure; hold on;
plot(ca, A, ca, Ap1, ca,Ap2, ca, Ap1p2);
legend('A', 'Ap1', 'Ap2', 'Ap1p2', 'Location', 'Best');
xlabel('Calcium'); ylabel('GLuR1 State Fractions');
title('Phosophorylation (NE applied)');
xlim([0,30]);

%recalculate ampa conducatance
gampa2 = A + 2*(Ap1+Ap2) + 4*(Ap1p2);

%plot augmented ampa conductance 
figure; hold on;
plot(ca, gampa1, ca, gampa2);
legend('No NE', 'NE applied', 'Location', 'Best');
xlabel('Calcium'); ylabel('AMPA Conductance');
xlim([0,30]);

% AMPA conductance vs frequency

figure
plot(f1,gampa2);
xlabel('Frequency [Hz]'); ylabel('Conductance');
title('AMPA Conductance vs. Frequency');


% AMPA Conductance with and without NE

figure
plot(f1,gampa1/gampa1(1),f1,gampa2./gampa1(1));
legend('No NE', 'With NE');
xlabel('Frequency [Hz]'); ylabel('Conductance');
title('AMPA Conductance vs. Frequency (with and without NE)');



%% Individual states

% Current
ErCa = -130         % [mV] reversal potential for Ca2+
I = 130.*(ge./gm);      % ge chosen such that ge/gm = 10^(-10)
                        % ge for different states chosen based on
                        % calculated calculated conductances above
                        
                        

% LTP

%Functions for phosphotase and kinase rates
EK = @(ca) 1+100.*(ca.^2)./(64 + (ca.^2));
EP = @(ca) 1+30.*(ca.^2)./(1+(ca.^2));

EP1 = EP(CaLTP);
EP2 = EP1;
EK1 = EK(CaLTP);
EK2 = EK1;

%Calculate and plot fraction in each state
A = (EP1.*EP2)./((EK2+EP2).*(EK1+EP1));
Ap1 = (EK1.*EP2)./((EK2+EP2).*(EK1+EP1));
Ap2 = (EK2.*EP1)./((EK2+EP2).*(EK1+EP1));
Ap1p2 = (EK1.*EK2)./((EK2+EP2).*(EK1+EP1));

%Model ampa conducatance (matches experimental)
gampaLTP = (A + 2*(Ap1+Ap2) + 4*(Ap1p2))./gampa1(1);

