% This script is an example for how to solve question 2 in the module for
% short-term synaptic plasticity of the G-node course.
% The differnce of this code with the file 'Section_2.m' is that here the
% built-in function fminsearch is used to find the parameters' values. 
% Written by Alex Loebel, 15.2.2014

%%
clear all;
% loading the data
load Data_depressing

% creating the mean trace
V=mean(V_depressing);
% dt=0.25;  % the sampling rate

%%
% fitting the synaptic and membrane time constants from the last EPSP. 
% first step is to isolate the last EPSP
y=V(4001:4200);
y=y-y(1);
t=0:size(y,2)-1;

% initial points of the search
tauin_est=3;
tau_est=30;
A_est=10;
x0 = [tauin_est,tau_est,A_est];

% calling the fminsearch function, which in turn uses the fit_epsp_func
% that I created, in which the difference between the estimated trace and
% experimental trace is calculated. The goal of the fminsearch function is
% therefore to find the best set of parameters that minimize this difference.  
[x,fval] = fminsearch(@(x) fit_epsp_func(x,y,t,dt),x0)

% the output values of fminsearch
tauinfit=x(1);
taufit=x(2);
A_optimal=x(3);

% examining how good is the fit
figure
plot(y)
hold on
plot(A_optimal.*(tauinfit/(tauinfit-taufit)).*(exp(-(t.*dt)./tauinfit)-exp(-(t.*dt)./taufit)),'k');
               %% v=A.*(tauin/(tauin-tau)).*(exp(-t./tauin*dt)-exp(-t./tau*dt)); %%how do i get


pause

%%
% fitting the parameters of the short-term depression
%%2i
% calculating the amplitudes from the experimental trace
% finding V0 and Vmax for the first eight EPSPs
for k=1:8
    V0_experiment(k)=V(400+200*(k-1));
    v_peak_experiment(k)=max(V(400+200*(k-1):400+200*k));
end

% finding V0 and Vmax for the last EPSP
V0_experiment(9)=V(4000);
v_peak_experiment(9)=max(V(4000:4200));
vdifferenceexperiment=v_peak_experiment-V0_experiment;

% initial conditions 
taurec_est=200;
U_est=0.2;
A_est=10;
x0 = [taurec_est, U_est, A_est];

% calling the fminsearch function, which in turn call the function I
% created to look for the differnce between the experimental and simulated
% amplitudes of the EPSPs. 
[x,fval] = fminsearch(@(x) fit_dynm_param_func(x,vdifferenceexperiment, tauinfit, taufit),x0)

% the output of the fminsearch function
taurecfit=x(1);
Ufit=x(2);
Afit=x(3);

%%
% calculating the amplitudes from the model's equations, so I can compare
% them to the experimental amplitudes and check how good the fit is. 
x_fit(1)=1;
for n=1:7
    x_fit(n+1)=x_fit(n)*(1-Ufit)*exp(-50/taurecfit)+1-exp(-50/taurecfit); %%3
end
x_fit(9)=x_fit(8)*(1-Ufit)*exp(-550/taurecfit)+1-exp(-550/taurecfit); %y -550 vs -50 = delta_t?
alpha_fit=Afit.*Ufit.*x_fit; %alpha parameter for easier writing/reading see 4

V0_fit(1)=0;  %%4
for k=1:7
    v_max_fit(k)=alpha_fit(k)*(alpha_fit(k)*taufit/(alpha_fit(k)*tauinfit-V0_fit(k)*(tauinfit-taufit)))^(taufit/(tauinfit-taufit)); %% 5
    V0_fit(k+1)=V0_fit(k).*exp(-50./taufit)+alpha_fit(k)*(tauinfit./(tauinfit-taufit)).*(exp(-50./tauinfit)-exp(-50./taufit));  %%6
end

v_max_fit(8)=alpha_fit(8)*(alpha_fit(8)*taufit/(alpha_fit(8)*tauinfit-V0_fit(8)*(tauinfit-taufit)))^(taufit/(tauinfit-taufit)); %%7
V0_fit(9)=V0_fit(8).*exp(-550./taufit)+alpha_fit(8)*(tauinfit./(tauinfit-taufit)).*(exp(-550./tauinfit)-exp(-550./taufit)); %%8
v_max_fit(9)=alpha_fit(9)*(alpha_fit(9)*taufit/(alpha_fit(9)*tauinfit-V0_fit(9)*(tauinfit-taufit)))^(taufit/(tauinfit-taufit)); %%9
vdifferencefit=v_max_fit-V0_fit; %%10

% plotting the results. 
figure
plot(vdifferenceexperiment,'k.-')
hold on
plot(vdifferencefit,'r.-')



