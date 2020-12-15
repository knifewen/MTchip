clear all;
close all;
clc;

N=160000;                     % number of points
fs=16e10;                   % sampling frequency
ts=1/fs;                    % sampling time
% t=(0:N-1)./fs;              % time
f=-fs/2:fs/N:fs/2-fs/N;     % operational frequency (Hz)
f_GHz = f./1e9;
interval = fs/N;

fc = 193e12;      %% Laser frequency

f1 = 5e9;        %% RF frequency
f2 = 5.01e9;
w_rf_1 = 2*pi*f1;

n=1.72;                      % n is Refractive index (double-stripe ng_si3n4=1.72, n_si2n4=1.9)
L_r=0.6966944;                   % L is Ring length in cm n1.72
L_mzi = 0.5.*L_r;

alpha=0.15;                 % alpha is Ring loss in dB/cm 
c=3e10;                     % c is speed of light in cm/s

psi1=1.*pi;                  % optical phase response (shift) 1  %%(positive shift to lower sideband (left))
psi2=pi;                  % optical phase response (shift) 2  already have the minus in () so the value of fshift is the position of ring 
psi3=pi;                  % optical phase response (shift) 2
psi_mzi=0.*pi; 

q1=exp(-1j*psi1);                % phase shift
q2=exp(-1j*psi2);                % phase shift
q3=exp(-1j*psi3);                % phase shift
q_mzi = exp(-1j*psi_mzi);                % phase shift


a_r = 10^(-alpha.*L_r./20);   %% round trip loss
a_wg = 10^(-alpha.*L_mzi./20);

% k_r1 = 0.96;  %%critical 0.012
% k_r2 = 0.86;
% k_r3 = 0.25;
k_r1 = 0.96;  %%critical 0.012
k_r2 = 0.68;
k_r3 = 0.25;

cc_r1 = sqrt(k_r1);
cc_r2 = sqrt(k_r2);
cc_r3 = sqrt(k_r3);

sc_r1 = sqrt(1-k_r1);
sc_r2 = sqrt(1-k_r2);
sc_r3 = sqrt(1-k_r3);

%% in freq
z_r = exp(1j.*2.*pi.*f.*n.*L_r./c);
z_mzi = exp(1j.*2.*pi.*f.*n.*L_mzi./c);
% %% in lambda
% z_r = exp(1j.*2.*pi.*n.*L_r./lambda);
% z_mzi = exp(1j.*2.*pi.*n.*L_mzi./lambda);

H_r1 = (sc_r1-a_r.*q1.*z_r.^(-1))./(1-sc_r1.*a_r.*q1.*z_r.^(-1));
H_r2 = (sc_r2-a_r.*q2.*z_r.^(-1))./(1-sc_r2.*a_r.*q2.*z_r.^(-1));
H_r3 = (sc_r3-a_r.*q3.*z_r.^(-1))./(1-sc_r3.*a_r.*q3.*z_r.^(-1));

H_wg = a_wg.*z_mzi.^(-1).*q_mzi;

H_r1_E = abs(H_r1);
H_r1_db = 20.*log10(abs(H_r1));
ph_r1 = -unwrap(angle(H_r1));


%% de-interleaver
%%% sc1, sc2 are self-coupling of first tunable coupler
%%% cc1, cc2 are cross-coupling of first tunable coupler
%%% phi1 is the phase shift of the first tunable coupler
%%% sc3, sc4 are self-coupling of second tunable coupler
%%% cc3, cc4 are cross-coupling of second tunable coupler
%%% phi2 is the phase shift of the second tunable coupler

k1 = 0.5;
k2 = 0.5;
k3 = 0.5;
k4 = 0.5;

cc1 = sqrt(k1);
cc2 = sqrt(k2);
cc3 = sqrt(k3);
cc4 = sqrt(k4);

sc1 = sqrt(1-k1);
sc2 = sqrt(1-k2);
sc3 = sqrt(1-k3);
sc4 = sqrt(1-k4);

phi1 = 0.5.*pi;
phi2 = 0.5.*pi;

H_bar = (sc3.*sc4.*exp(-1j.*phi2)-cc3.*cc4).*(sc1.*sc2.*exp(-1j.*phi1)-cc1.*cc2).*H_r1.*H_r3...
        -(sc4.*cc3.*exp(-1j.*phi2)-sc3.*cc4).*(sc1.*cc2.*exp(-1j.*phi1)-sc2.*cc1).*H_r2.*H_wg;
    
H_cross = (sc3.*sc4.*exp(-1j.*phi2)-cc3.*sc4).*(cc1.*sc2.*exp(-1j.*phi1)-cc1.*cc2).*H_r1.*H_r3...
        -(sc3.*sc4-cc3.*cc4.*exp(-1j.*phi2)).*(sc1.*cc2.*exp(-1j.*phi1)-sc2.*cc1).*H_r2.*H_wg;

%%% amplitude
H_bar_db = 20.*log10(abs(H_bar));
H_cross_db = 20.*log10(abs(H_cross));
%%% phase
ph_bar = -unwrap(angle(H_bar));
ph_cross = -unwrap(angle(H_cross));

%% UMZI
H_MZI_bar = (sc3.*sc4.*exp(-1j.*phi2)-cc3.*cc4).*(sc1.*sc2.*exp(-1j.*phi1)-cc1.*cc2)...
        -(sc4.*cc3.*exp(-1j.*phi2)-sc3.*cc4).*(sc1.*cc2.*exp(-1j.*phi1)-sc2.*cc1).*H_wg;
    
H_MZI_cross = (sc3.*sc4.*exp(-1j.*phi2)-cc3.*sc4).*(cc1.*sc2.*exp(-1j.*phi1)-cc1.*cc2)...
        -(sc3.*sc4-cc3.*cc4.*exp(-1j.*phi2)).*(sc1.*cc2.*exp(-1j.*phi1)-sc2.*cc1).*H_wg;

%%% amplitude
H_MZI_bar_db = 20.*log10(abs(H_MZI_bar));
H_MZI_cross_db = 20.*log10(abs(H_MZI_cross));
%%% phase
ph_MZI_bar = -unwrap(angle(H_MZI_bar));
ph_MZI_cross = -unwrap(angle(H_MZI_cross));
%% plot ring
figure;

plot(f_GHz,H_r1_db,'linewidth',2);
ylabel('Transmission (dB)')
grid on
hold on
ax=gca;
ax.FontSize=16;
ay=gca;
ay.FontSize=16;
% xlim([-15 15]);
% title('wave shaper');
% xlabel('detuning phi');
xlabel('Frequency (GHz)');
%legend('OCB','+2OSB','-2OSB');
% %% plot ring
% figure;
% 
% plot(f_GHz,ph_r1./pi,'linewidth',2);
% ylabel('Phase (\pi)')
% grid on
% hold on
% ax=gca;
% ax.FontSize=16;
% ay=gca;
% ay.FontSize=16;
% % xlim([-15 15]);
% % title('wave shaper');
% % xlabel('detuning phi');
% xlabel('Frequency (GHz)');
% %legend('OCB','+2OSB','-2OSB');
%% plot interleaver
%% plot bar
figure;

plot(f_GHz,H_bar_db,'linewidth',2);
hold on;
plot(f_GHz,H_cross_db,'linewidth',2);
ylabel('Transmission (dB)')
grid on
hold on
ax=gca;
ax.FontSize=16;
ay=gca;
ay.FontSize=16;
% xlim([-15 15]);
% title('wave shaper');
% xlabel('detuning phi');
xlabel('Frequency (GHz)');
legend('bar','cross');

%%% phase
figure;
plot(f_GHz,ph_bar,'linewidth',2);
hold on
plot(f_GHz,ph_cross,'linewidth',2);
ylabel('Phase (\pi)')
grid on
hold on
ax=gca;
ax.FontSize=16;
ay=gca;
ay.FontSize=16;
% xlim([-15 15]);
% title('wave shaper');
% xlabel('detuning phi');
xlabel('Frequency (GHz)');
legend('bar','cross');

% %% plot cross
% figure;
% 
% plot(f_GHz,H_cross_db,'linewidth',2);
% ylabel('Transmission (dB)')
% grid on
% hold on
% ax=gca;
% ax.FontSize=16;
% ay=gca;
% ay.FontSize=16;
% % xlim([-15 15]);
% % title('wave shaper');
% % xlabel('detuning phi');
% xlabel('Frequency (GHz)');
% %legend('OCB','+2OSB','-2OSB');
% 
% %%%phase
% figure;
% plot(f_GHz,ph_cross,'linewidth',2);
% ylabel('Phase (\pi)')
% grid on
% hold on
% ax=gca;
% ax.FontSize=16;
% ay=gca;
% ay.FontSize=16;
% % xlim([-15 15]);
% % title('wave shaper');
% % xlabel('detuning phi');
% xlabel('Frequency (GHz)');
% %legend('OCB','+2OSB','-2OSB');
% figure
% plot(abs(H_bar));
% 
% %% plot MZI cross
% figure;
% 
% plot(f_GHz,H_MZI_cross_db,'linewidth',2);
% ylabel('Transmission (dB)')
% grid on
% hold on
% ax=gca;
% ax.FontSize=16;
% ay=gca;
% ay.FontSize=16;
% % xlim([-15 15]);
% % title('wave shaper');
% % xlabel('detuning phi');
% xlabel('Frequency (GHz)');
% %legend('OCB','+2OSB','-2OSB');
% 
% %%%phase
% figure;
% plot(f_GHz,ph_MZI_cross,'linewidth',2);
% ylabel('Phase (\pi)')
% grid on
% hold on
% ax=gca;
% ax.FontSize=16;
% ay=gca;
% ay.FontSize=16;
% % xlim([-15 15]);
% % title('wave shaper');
% % xlabel('detuning phi');
% xlabel('Frequency (GHz)');
% %legend('OCB','+2OSB','-2OSB');

