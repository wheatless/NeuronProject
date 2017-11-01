% 
clear all
close all
clc

% Constants
% Defined according to Shouval 2001
f = 100;            % [Hz] Stimulation frequency
tauCa = 1;          % Decay time constant of Ca2+ in spine
tauf = 1;           % Fast: Time constant of Ca2+ current thru NMDAR
taus = 1;           % Slow: Time constant of Ca2+ current thru NMDAR
Mg = 1;             % [mM] Magnesium concentration 
Vr = 130;           % [mV] Reversal potential of Ca2+
a1 = (100-0)/(30-(-70));         % 1 [Hz/mV] Proportionality constant relating stim frequency to voltage
b = -70;                         % [mV] Y-intercept of V(f) curve

V = @(f) a1.*f + b;              % [mV] Postsynaptic potential
%Nf = func of fast NMDA rec;             % Fast: Magnitude of Ca2+ current thru NMDAR
%Ns = func of slow NMDA rec;             % Slow: Magnitude of Ca2+ current thru NMDAR
Nf = 1; Ns = 1;

v = V(f);

B = @(V) 1./(1 + exp(-0.062*V).*(Mg/3.57));
H = @(V) B(V).*(V-Vr);
Gnmda = tauCa.*(tauf.*Nf + taus.*Ns);   % Gain of Ca2+ influx thru NMDA rec. They use 
                                        % Gnmda = 0.01 and 0.03

CaSS = H(v).*f.*Gnmda;
%%
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
