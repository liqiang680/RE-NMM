%% Function to generate spike discharges in Rolandic epilepsy (RE)
% This code is about the neural mass model RE-NMM. 
% It is constructed to study the competitive dynamics observed between spikes and spindles in Rolandic epilepsy (RE).
% NB: for different states the Mg2+ block factor, conductance of h-type current, inhibitory projection strengths from TRN, and conductance of T-type current must be changed.

% Reference: code Costa2016_model
% M. S. Costa and A. Weigenand, et al. 
% "A thalamocortical neural mass model of the eeg during nrem sleep and its response to auditory stimulation." 
% PLOS Computational Biology, vol.12, no.9, 2016.

%%
% 
function Rolandic_model()

clear all; 
clc;
% Parameters
% Membrane capacitance (uF/cm2)
Cm=1;
% Membrane time constant (ms)
taop=30; 
taoi=30; 
taot=20; 
taor=20;
% Maximal firing rate (ms)
Qp_max=30e-3; 
Qi_max=60e-3; 
Qt_max=400e-3; 
Qr_max=400e-3;
% Firing threshold (mV)
theta_p=-58.5; 
theta_i=-58.5; 
theta_t=-58.5; 
theta_r=-58.5;
% Inverse neural gain (mV)
sigma_p=6;
sigma_i=6; 
sigma_t=6; 
sigma_r=6;
% Synaptic rate constant (ms-1)
gammae=70e-3;
gammag=100e-3;  
gamma_n=30e-3;
gamma_N=100e-3;    
gammagp=58.6e-3;

% Mg2+ block factor (mM)
mu=2.2; % reflecting the density of NMDA current

% Connectivity constant
Nip=0; 
Npi=0;
Nii=0;
Nit=0; 
Npp=70;
Nt1p=3;
Nr1p=2.6;
Nr1t1=4;

% The inhibitory projections from TRN
Nt1r1=5.5; % from PV to VPM
Nt2r1=2; % from PV to VPL
Nt2r2=1.9; % from SOM to VPL

Nrr=15;
Npt1=3; 
Nrr2=25; 
Nr2t2=1.5; 
Nr1r2=2.1; 
Npt2=2; 
Nt2p=2.1; 
Nr1t2=3; 
Nr2r1=0;  
Nt2t1=1; 
Nt1t2=1;
% Input rate of synaptic channel (ms)
WAMPA=1; 
WGABA=1;
WNMDA=1;

% Conductivity of ion channel
g_T_t=2.6;  % the conductance of T-type Calcium current (VPM/VPL)
g_T_r1=1.6; %  (PV)
g_T_r2=1.6; % (SOM)

g_L_K=0.042; 
g_L_Kr2=0.042; 
g_L_Kr1=0.042;

% Conductance of h-type current
g_h=0.056;  % reflecting the density of h-type current

% Nernst reversal potential (mV)
E_L_p=-64;
E_L_i=-64;
E_L_t=-70;
E_L_r=-70; 
E_K=-100;
E_Ca_t=120; 
E_Ca_r1=120; 
E_Ca_r2=120; 
EAMPA=0; 
EGABA=-70; 
E_h=-40; 
% Standard deviation of cortical/thalamic background noise (ms-1)
phi_C_sd=120e-3;
phi_T_sd=20e-3;
phi_S_sd=10e-3;
% Calcium dynamic
alphaCa=-51.8e-6;tao_Ca=10;Ca0=2.4e-4;
% h-type current
k1=2.5e7;k2=4e-4;k3=1e-1;k4=1e-3;np=4;ginc=2;
% axonal rate
v=120e-3;

% Time (ms-1)
start = 0.0;
stop = 61500 ; 
dt = 0.1;
time_array = [start:dt:stop];
vec_len = length(time_array);

% Outputs
Vp = zeros(1,vec_len);
Vi = zeros(1,vec_len);
Vt = zeros(1,vec_len);
Vr = zeros(1,vec_len);
Sep = zeros(1,vec_len);
Sepdot = zeros(1,vec_len);
Sei = zeros(1,vec_len);
Seidot = zeros(1,vec_len);
Set = zeros(1,vec_len);
Setdot = zeros(1,vec_len);
Ser = zeros(1,vec_len);
Serdot = zeros(1,vec_len);
Serplus = zeros(1,vec_len);
Serplusdot = zeros(1,vec_len);
Sgp = zeros(1,vec_len);
Sgpdot = zeros(1,vec_len);
Sgi = zeros(1,vec_len);
Sgidot = zeros(1,vec_len);
Srt = zeros(1,vec_len);
Srtdot = zeros(1,vec_len);
Srr = zeros(1,vec_len);
Srrdot = zeros(1,vec_len);
phip = zeros(1,vec_len);
phipdot = zeros(1,vec_len);
phit1 = zeros(1,vec_len);
phitdot1 = zeros(1,vec_len);
phit2 = zeros(1,vec_len);
phitdot2 = zeros(1,vec_len);
h_T_t = zeros(1,vec_len);
h_T_t2 = zeros(1,vec_len);
h_T_r = zeros(1,vec_len);
h_T_r2 = zeros(1,vec_len);
m_h1 = zeros(1,vec_len);
m_h2 = zeros(1,vec_len);
m_h12 = zeros(1,vec_len);
m_h22 = zeros(1,vec_len);
concentration_Ca = zeros(1,vec_len);
concentration_Ca2 = zeros(1,vec_len);
Vr2= zeros(1,vec_len);
Ser2= zeros(1,vec_len);
Serdot2= zeros(1,vec_len);
Srr2= zeros(1,vec_len);
Srrdot2= zeros(1,vec_len);
Vt2=zeros(1,vec_len);
Set2 = zeros(1,vec_len);
Setdot2= zeros(1,vec_len);
Srt2= zeros(1,vec_len);
Srtdot2 = zeros(1,vec_len);
Srt2_slow= zeros(1,vec_len);
Srt2_slowdot2= zeros(1,vec_len);
Srr_slow= zeros(1,vec_len);
Srr_slowdot= zeros(1,vec_len);
Ser2plus= zeros(1,vec_len);
Ser2plusdot2= zeros(1,vec_len);

% Initialize outputs
Ser2plus(1) = 0;
Ser2plusdot2(1) = 0;
Srt2_slow(1) = 0;
Srt2_slowdot2(1) = 0;
Srr_slow(1) = 0;
Srr_slowdot(1) = 0;
Vp(1) = -64;
Vi(1) = -64;
Vt(1) = -68;
Vr(1) = -68;
Serplus(1) = 0;
Serplusdot(1) = 0;
Vt2(1) = -68;
Set2(1) = 0;
Setdot2(1) = 0;
Srt2(1) = 0;
Srtdot2(1) = 0;
Vr2(1) = -68;
Ser2(1) = 0;
Serdot2(1) = 0;
Srr2(1) = 0;
Srrdot2(1) = 0;
Sep(1) = 0;
Sepdot(1) = 0;
Sei(1) = 0;
Seidot(1) = 0;
Set(1) = 0;
Setdot(1) = 0;
Ser(1) = 0;
Serdot(1) = 0;
Sgp(1) = 0;
Sgpdot(1) = 0;
Sgi(1) = 0;
Sgidot(1) = 0;
Srt(1) = 0;
Srtdot(1) = 0;
Srr(1) = 0;
Srrdot(1) = 0;
phip(1) = 0;
phipdot(1) = 0;
phit1(1) = 0;
phitdot1(1) = 0;
phit2(1) = 0;
phitdot2(1) = 0;
h_T_t(1) = 0;
h_T_t2(1) = 0;
h_T_r(1) = 0;
h_T_r2(1) = 0;
m_h1(1) = 0;
m_h2(1) = 0;
m_h12(1) = 0;
m_h22(1) = 0;
concentration_Ca(1) = 2.4E-4;
concentration_Ca2(1) = 2.4E-4;

% Noise 
noise1 = gammae^2*sqrt(phi_C_sd*dt)*randn(1,vec_len); 
noise2 = gammae^2*sqrt(phi_C_sd*dt)*randn(1,vec_len); 
noise3 = gammae^2*sqrt(phi_T_sd*dt)*randn(1,vec_len); 
noise4 = gammae^2*sqrt(phi_S_sd*dt)*randn(1,vec_len); 

for t = 2:(stop/dt)

% Firing rate
C=(pi./sqrt(3));
Qp=Qp_max./(1+exp(-C*(Vp(t-1)-theta_p)/(sigma_p)));
Qi=Qi_max./(1+exp(-C*(Vi(t-1)-theta_i)/(sigma_i)));
Qt1=Qt_max./(1+exp(-C*(Vt(t-1)-theta_t)/(sigma_t)));
Qt2=Qt_max./(1+exp(-C*(Vt2(t-1)-theta_t)/(sigma_t)));
Qr1=Qr_max./(1+exp(-C*(Vr(t-1)-theta_r)/(sigma_r)));
Qr2=Qr_max./(1+exp(-C*(Vr2(t-1)-theta_r)/(sigma_r)));
 
% Leak current
J_L_p=(Vp(t-1)-E_L_p);
J_L_i=(Vi(t-1)-E_L_i);
J_L_t=(Vt(t-1)-E_L_t);
J_L_t2=(Vt2(t-1)-E_L_t);
J_L_r=(Vr(t-1)-E_L_r);
J_L_r2=(Vr2(t-1)-E_L_r);

% Synaptic current
JAMPA_Sep=WAMPA*Sep(t-1)*(Vp(t-1)-EAMPA);
JAMPA_Sei=WAMPA*Sei(t-1)*(Vi(t-1)-EAMPA);
JAMPA_Set=WAMPA*Set(t-1)*(Vt(t-1)-EAMPA);
JAMPA_Set2=WAMPA*Set2(t-1)*(Vt2(t-1)-EAMPA);
JAMPA_Ser=WAMPA*Ser(t-1)*(Vr(t-1)-EAMPA);
JAMPA_Ser2=WAMPA*Ser2(t-1)*(Vr2(t-1)-EAMPA);

JGABA_Sgp=WGABA*Sgp(t-1)*(Vp(t-1)-EGABA);
JGABA_Sgi=WGABA*Sgi(t-1)*(Vi(t-1)-EGABA);
JGABA_Srt=WGABA*Srt(t-1)*(Vt(t-1)-EGABA);
JGABA_Srt2=WGABA*(Srt2(t-1))*(Vt2(t-1)-EGABA);
JGABA_Srr=WGABA*(Srr(t-1))*(Vr(t-1)-EGABA);
JGABA_Srr2=WGABA*Srr2(t-1)*(Vr2(t-1)-EGABA);

fun=1./(1+exp(-mu.*0.062*Vr(t-1))*(1./3.57));
JNMDA_Serplus=WNMDA*fun*Serplus(t-1)*(Vr(t-1)-EAMPA); %NMDA current
fun2=1./(1+exp(-mu.*0.062*Vr2(t-1))*(1./3.57));
JNMDA_Ser2plus=WNMDA*fun2*Ser2plus(t-1)*(Vr2(t-1)-EAMPA); %NMDA current

% Activation function 
m_limit_t=1./(1+exp(-(Vt(t-1)+59)/6.2));
m_limit_t2=1./(1+exp(-(Vt2(t-1)+59)/6.2));
m_limit_r=1./(1+exp(-(Vr(t-1)+52)/7.4)); %  
m_limit_r2=1./(1+exp(-(Vr2(t-1)+55.5)/7.4));
m_limit_h=1./(1+exp((Vt(t-1)+75)/5.5));
m_limit_h2=1./(1+exp((Vt2(t-1)+75)/5.5));
% inactivation function 
h_limit_t=1./(1+exp((Vt(t-1)+81)/4));
h_limit_t2=1./(1+exp((Vt2(t-1)+81)/4));
h_limit_r=1./(1+exp((Vr(t-1)+80)/5));
h_limit_r2=1./(1+exp((Vr2(t-1)+80)/5));

%Intrinsic currents
I_LK_t=g_L_K*(Vt(t-1)-E_K);
I_LK_t2=g_L_K*(Vt2(t-1)-E_K);
I_LK_r=g_L_Kr1*(Vr(t-1)-E_K);
I_LK_r2=g_L_Kr2*(Vr2(t-1)-E_K);
I_T_t=g_T_t*m_limit_t*m_limit_t*h_T_t(t-1).*(Vt(t-1)-E_Ca_t);
I_T_t2=g_T_t*m_limit_t2*m_limit_t2*h_T_t2(t-1).*(Vt2(t-1)-E_Ca_t);
I_T_r=g_T_r1*m_limit_r*m_limit_r*h_T_r(t-1)*(Vr(t-1)-E_Ca_r1);
I_T_r2=g_T_r2*m_limit_r2*m_limit_r2*h_T_r2(t-1)*(Vr2(t-1)-E_Ca_r2);
tao_h_t2=(30.8+(211.4+exp((Vt2(t-1)+115.2)/5))./(1+exp((Vt2(t-1)+86)/3.2)))/3.7371928;
tao_h_r2=(85+1./(exp((Vr2(t-1)+48)/4)+exp(-(Vr2(t-1)+407)/50)))/3.7371928;
tao_h_t=(30.8+(211.4+exp((Vt(t-1)+115.2)/5))./(1+exp((Vt(t-1)+86)/3.2)))/3.7371928;
tao_h_r=(85+1./(exp((Vr(t-1)+48)/4)+exp(-(Vr(t-1)+407)/50)))/3.7371928;
h_T_t(t)=h_T_t(t-1)+dt*((h_limit_t-h_T_t(t-1))./tao_h_t);
h_T_t2(t)=h_T_t2(t-1)+dt*((h_limit_t2-h_T_t2(t-1))./tao_h_t2);


tao_m_h=(20+1000./(exp((Vt(t-1)+71.5)/14.2)+exp(-(Vt(t-1)+89)/11.6)));
tao_m_h2=(20+1000./(exp((Vt2(t-1)+71.5)/14.2)+exp(-(Vt2(t-1)+89)/11.6)));
Ph=k1*(concentration_Ca(t-1))^np./(k1*(concentration_Ca(t-1))^np+k2);
Ph2=k1*(concentration_Ca2(t-1))^np./(k1*(concentration_Ca2(t-1))^np+k2);
h_T_r(t)=h_T_r(t-1)+dt*((h_limit_r-h_T_r(t-1))./tao_h_r);
h_T_r2(t)=h_T_r2(t-1)+dt*((h_limit_r2-h_T_r2(t-1))./tao_h_r2);
m_h1(t)=m_h1(t-1)+dt*((m_limit_h*(1-m_h2(t-1))-m_h1(t-1))./tao_m_h-k3*Ph*m_h1(t-1)+k4*m_h2(t-1));
m_h2(t)=m_h2(t-1)+dt*(k3*Ph*m_h1(t-1)-k4*m_h2(t-1));
m_h12(t)=m_h12(t-1)+dt*((m_limit_h2*(1-m_h22(t-1))-m_h12(t-1))./tao_m_h2-k3*Ph2*m_h12(t-1)+k4*m_h22(t-1));
m_h22(t)=m_h22(t-1)+dt*(k3*Ph2*m_h12(t-1)-k4*m_h22(t-1));
concentration_Ca(t)=concentration_Ca(t-1)+dt*(alphaCa*I_T_t-(concentration_Ca(t-1)-Ca0)./tao_Ca);
concentration_Ca2(t)=concentration_Ca2(t-1)+dt*(alphaCa*I_T_t2-(concentration_Ca2(t-1)-Ca0)./tao_Ca);
I_h=g_h*(m_h1(t-1)+ginc*m_h2(t-1))*(Vt(t-1)-E_h); %h-type current
I_h2=g_h*(m_h12(t-1)+ginc*m_h22(t-1))*(Vt2(t-1)-E_h);  %h-type current


% Membrane potential 
Vp(t)=Vp(t-1)+dt*((-J_L_p-JAMPA_Sep-JGABA_Sgp)./taop);
Vi(t)=Vi(t-1)+dt*((-J_L_i-JAMPA_Sei-JGABA_Sgi)./taoi);
Vt(t)=Vt(t-1)+dt*((-J_L_t-JAMPA_Set-JGABA_Srt-Cm*taot*(I_LK_t+I_T_t+I_h))./taot);
Vt2(t)=Vt2(t-1)+dt*((-J_L_t2-JAMPA_Set2-JGABA_Srt2-Cm*taot*(I_LK_t2+I_T_t2+I_h2))./taot);
Vr(t)=Vr(t-1)+dt*((-J_L_r-JAMPA_Ser-JGABA_Srr-JNMDA_Serplus-Cm*taor*(I_LK_r+I_T_r))./taor);
Vr2(t)=Vr2(t-1)+dt*((-J_L_r2-JAMPA_Ser2-JGABA_Srr2-JNMDA_Ser2plus-Cm*taor*(I_LK_r2+I_T_r2))./taor);

% Proportion of open synaptic channels
Sep(t) = Sep(t-1) + Sepdot(t-1)*dt;
Sepdot(t) = Sepdot(t-1) + dt*(gammae^2*(Npp*Qp+Npt1*phit1(t-1)+Npt2*phit2(t-1)-Sep(t-1))-2*gammae*Sepdot(t-1) )+noise1(t-1);

Sei(t) = Sei(t-1) + Seidot(t-1)*dt;
Seidot(t) = Seidot(t-1) + dt*(gammae^2*(Nip*Qp+Nit*1-Sei(t-1))-2*gammae*Seidot(t-1) )+noise2(t-1);

Set(t) = Set(t-1) + Setdot(t-1)*dt;
Setdot(t) = Setdot(t-1) + dt*(gammae^2*(Nt1p*phip(t-1)+Nt2t1*Qt2-Set(t-1))-2*gammae*Setdot(t-1) )+noise3(t-1);

Set2(t) = Set2(t-1) + Setdot2(t-1)*dt;
Setdot2(t) = Setdot2(t-1) + dt*(gammae^2*(Nt2p*phip(t-1)+Nt1t2*Qt1-Set2(t-1))-2*gammae*Setdot2(t-1));

Ser(t) = Ser(t-1) + Serdot(t-1)*dt;
Serdot(t) = Serdot(t-1) + dt*(gammae^2*(Nr1t1*Qt1+Nr1t2*Qt2+Nr1p*phip(t-1)-Ser(t-1))-2*gammae*Serdot(t-1) );

Serplus(t) = Serplus(t-1) + Serplusdot(t-1)*dt;
Serplusdot(t) = Serplusdot(t-1) + dt*(gammag*gamma_n*(Nr1t1*Qt1+Nr1t2*Qt2+Nr1p*phip(t-1)-Serplus(t-1))-2*gamma_n*Serplusdot(t-1) );

Ser2(t) = Ser2(t-1) + Serdot2(t-1)*dt;
Serdot2(t) = Serdot2(t-1) + dt*(gammae.*gammae.*(Nr2t2*Qt2-Ser2(t-1))-2*gammae*Serdot2(t-1) )+noise4(t-1) ;

Ser2plus(t) = Ser2plus(t-1) + Ser2plusdot2(t-1)*dt;
Ser2plusdot2(t) = Ser2plusdot2(t-1) + dt*(gammag*gamma_n*(Nr2t2*Qt2-Ser2plus(t-1))-2*gamma_n*Ser2plusdot2(t-1));

Sgp(t) = Sgp(t-1) + Sgpdot(t-1)*dt;
Sgpdot(t) = Sgpdot(t-1) + dt*(gammagp^2*(Npi*Qi-Sgp(t-1))-2*gammagp*Sgpdot(t-1) );

Sgi(t) = Sgi(t-1) + Sgidot(t-1)*dt;
Sgidot(t) = Sgidot(t-1) + dt*(gammagp^2*(Nii*Qi-Sgi(t-1))-2*gammagp*Sgidot(t-1) );

Srt(t) = Srt(t-1) +Srtdot(t-1)*dt;
Srtdot(t) = Srtdot(t-1) + dt*(gamma_N.*gamma_N*(Nt1r1*Qr1-Srt(t-1))-2*gamma_N*Srtdot(t-1) );

Srt2(t) = Srt2(t-1) +Srtdot2(t-1)*dt;
Srtdot2(t) = Srtdot2(t-1) + dt*(gamma_N.*gamma_N*(Nt2r1*Qr1-Srt2(t-1)) +gamma_N.*gamma_N*Nt2r2*Qr2  -2*gamma_N*Srtdot2(t-1));

Srt2_slow(t) = Srt2_slow(t-1) +Srt2_slowdot2(t-1)*dt;
Srt2_slowdot2(t) = Srt2_slowdot2(t-1) + dt*(gamma_N^2*(Nt2r2*Qr2-Srt2_slow(t-1))-2*gamma_N*Srt2_slowdot2(t-1) );

Srr(t) = Srr(t-1) +Srrdot(t-1)*dt;
Srrdot(t) = Srrdot(t-1) + dt*(gamma_N^2*(Nrr*Qr1-Srr(t-1)+Nr1r2*Qr2)-2*gamma_N*Srrdot(t-1) );

Srr_slow(t) = Srr_slow(t-1) +Srr_slowdot(t-1)*dt;
Srr_slowdot(t) = Srr_slowdot(t-1) + dt*(gamma_N^2*(Nr1r2*Qr2-Srr_slow(t-1))-2*gamma_N*Srr_slowdot(t-1) );

Srr2(t) = Srr2(t-1) +Srrdot2(t-1)*dt;
Srrdot2(t) = Srrdot2(t-1) + dt*(gamma_N.*gamma_N*(Nrr2*Qr2+Nr2r1*Qr1-Srr2(t-1))-2*gamma_N*Srrdot2(t-1) );

% Delayed firing rate
phip(t) = phip(t-1) +phipdot(t-1)*dt;
phipdot(t) = phipdot(t-1) + dt*(v^2*(Qp-phip(t-1))-2*v*phipdot(t-1));

phit1(t) = phit1(t-1) +phitdot1(t-1)*dt;
phitdot1(t) = phitdot1(t-1) + dt*(v^2*(Qt1-phit1(t-1))-2*v*phitdot1(t-1));

phit2(t) = phit2(t-1) +phitdot2(t-1)*dt;
phitdot2(t) = phitdot2(t-1) + dt*(v^2*(Qt2-phit2(t-1))-2*v*phitdot2(t-1));
end

time_array2=start:dt/1000:stop/1000;
V1=Vp(15000-1:615000);
figure;
plot(time_array2(15000-1:615000),V1),title('Vp');
xlim([1.5,61.5]);
