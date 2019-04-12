%% Stephen Lisius
%% Parameter input for ablate function

clear
close all

%inputs
%laser properties
fx=1.7e12;     %laser power flux (w/m^2) (J/m^2 s)
pl=6e-9;     %pulse duration (s)
pf=20;          %pulse frequency (hz)
ot=1e-6;    %off time between pulses
bsd=.1;    %beam spot diameter (m)
N=3000;        %number of nodes
dx=1e-9;         %thickness of each layer (m)
dt_ab=1e-14;            %time step (s)
dt_cd=1e-9;
Ti=300;         %initial temperature
dmass=1;        %debris mass (kg)

x=ones(N,3);
x(:,1)=Ti;
x(:,2)=0;
x(:,3)=2;
x(1,3)=0;
x(1,1)=0;

[imp_pp, x, normmass]=ablate(fx, pl, bsd, N, dx, dt_ab, x);

figure
plot(x(:,1))

dv_pp=imp_pp/dmass;
dv_ps=dv_pp*pf;


% [x, tstore2]=cooldown(ot, N, dx, dt_cd, bsd, x);
% 
% figure
% plot(x(:,1))