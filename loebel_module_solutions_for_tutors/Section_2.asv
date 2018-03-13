% This script is an example for how to solve question 2 in the module for
% short-term synaptic plasticity of the G-node course.
% Written by Alex Loebel, 15.2.2014
clear all;
close all;
% loading the data to analyse
load Data_depressing

% creating the mean trace, as we analyse the only the average of the
% response
V= mean(V_depressing);

figure
plot(V,'.-')

dt=0.25;  % the sampling rate of 4kHz leads to time difference of 0.25ms
% between measurements points


% calculating the memebrane time constant and the synaptic inactivation
% time constant.
% The following is a self written minimazation algorithm for the error.
% There are many other ways of doing this, e.g., using built in functions
% like fminsearch.
y=V(4001:4120);  % using the ninth EPSP for the analysis
y=y-y(1);  % so the initial point of the EPSP that we analyse will be zero.

figure
plot(y,'.-')

t=0:size(y,2)-1;  % 'fake' time, later it is multiplied by dt to get the real
% value of the time.
c=10^9;    % initial error value. to minimize later.

tau_optimal=0;
tauin_optimal=0;
A_optimal=0;

b=0;

num_of_intervals=40; % the size of the grid
% range of values for the different parameters.
tau_min=1;
tau_max=50;
dtau=(tau_max-tau_min)/num_of_intervals;  % initial grid step for the membrane
% time constant.

tauin_min=1;
tauin_max=8;
dtauin=(tauin_max-tauin_min)/num_of_intervals;  % initial grid step for the inactivation
% time constant.

A_min=1;
A_max=80;
dA=(A_max-A_min)/num_of_intervals; % initial grid step for the amplitude dummy variable.

kk=1;
while kk<=4
    % the following flags are used for boundry conditions of the different
    % variables, e.g., if we get that tauin is larger than the max value of
    % 8ms, the max value is set to a higher value. 
    aa=0;
    aa1=0;
    aa2=0;
    aa3=0;
    aa4=0;
    aa5=0;
    
    for tau=[tau_min:dtau:tau_max]
        for tauin=[tauin_min:dtauin:tauin_max]
            for A=[A_min:dA:A_max]
                % calculating the simulated trace
                v=A.*(tauin/(tauin-tau)).*(exp(-t./tauin*dt)-exp(-t./tau*dt));
                
                % calculating the least-mean-squares between the
                % experimental and the simulated traces. 
                b=sum((y-v).^2);
                if b<c
                    c=b;
                    tau_optimal=tau;
                    tauin_optimal=tauin;
                    A_optimal=A;
                end
            end
        end
    end

    % determining the new borders of the the different parameters.
    % at each iteration of the fitting, the resolution of the grid becomes
    % finer and finer.
    if  tau_optimal==tau_min
        tau_min=max(0,tau_min-2*dtau);
        tau_max=tau_max-2*dtau;
        aa=1;
    end
    if tau_optimal==tau_max
        tau_max=tau_max+2*dtau;
        tau_min=tau_min+2*dtau;
        aa1=1;
    end
    if tauin_optimal==tauin_min
        tauin_min=max(0,tauin_min-2*dtauin);
        tauin_max=tauin_max-2*dtauin;
        aa2=1;
    end
    if tauin_optimal==tauin_max
        tauin_max=tauin_max+2*dtauin;
        tauin_min=tauin_min+2*dtauin;
        aa3=1;
    end
    if A_optimal==A_min
        A_min=max(0,A_min-2*dA);
        A_max=A_max-2*dA;
        aa4=1;
    end
    if A_optimal==A_max
        A_max=A_max+2*dA;
        A_min=A_min+2*dA;
        aa5=1;
    end

    % making sure the boundries are above 0 for each of the
    % parameters.
    if (aa==0 & aa1==0 & aa2==0 & aa3==0 & aa4==0 & aa5==0)
        tau_min=tau_optimal-2*dtau;
        tau_max=tau_optimal+2*dtau;
        dtau=(tau_max-tau_min)/20;
        tauin_min=tauin_optimal-2*dtauin;
        tauin_max=tauin_optimal+2*dtauin;
        dtauin=(tauin_max-tauin_min)/20;
        A_min=A_optimal-2*dA;
        A_max=A_optimal+2*dA;
        dA=(A_max-A_min)/20;
        kk=kk+1
    end
end

% plotting the fit, to check how good it is. A button should be pressed to
% release Matlab from the pause mode. 
figure
plot(y)
hold on
plot(A_optimal.*(tauin_optimal/(tauin_optimal-tau_optimal)).*(exp(-t./tauin_optimal*dt)-exp(-t./tau_optimal*dt)),'k');
A_optimal
pause

% finding V0 and Vmax for the first eight EPSPs
for k=1:8
    V0_experiment(k)=V(400+200*(k-1));
    v_peak_experiment(k)=max(V(400+200*(k-1):400+200*k));
end

% finding V0 and Vmax for the last EPSP
V0_experiment(9)=V(4000);
v_peak_experiment(9)=max(V(4000:4200));
vdifferenceexperiment=v_peak_experiment-V0_experiment;

% considering the values for the membrane time constant and the
% inactivation time constant that were found above
tauin=tauin_optimal;
tau=tau_optimal;

% the following is the same fitting process as used above for finding 
% the membrane time constant and the inactivation time constant, only here
% it is adapted for the dynamic equations and for finding the three
% parameters of utilization, recovery time constant and the absolute
% efficacy. 

% the error we want to minimize
c=10^9;
b=0;

taurec_optimal=0;
U_optimal=0;
A_optimal=0;

taurec_min=250;
taurec_max=1500;
dtaurec=(taurec_max-taurec_min)/num_of_intervals;

U_min=0.1;
U_max=1;
dU=(U_max-U_min)/num_of_intervals;

A_min=10;
A_max=100;
dA=(A_max-A_min)/num_of_intervals;

kk=1
while kk<=4
    ee=0;
    ee1=0;
    bb=0;
    bb1=0;
    cc=0;
    dd=0;
    
    % fitting
    for U_fit=[U_min:dU:U_max]
        for A_fit=[A_min:dA:A_max]
            for taurec_fit=[taurec_min:dtaurec:taurec_max]
%%delta T1 = 50 delta T2 = 550
                x_fit(1)=1;
                for n=1:8
                    x_fit(n+1)=x_fit(n)*(1-U_fit)*exp(-50/taurec_fit)+1-exp(-50/taurec_fit);
                end
                x_fit(9)=x_fit(8)*(1-U_fit)*exp(-550/taurec_fit)+1-exp(-550/taurec_fit);
                alpha_fit=A_fit.*U_fit.*x_fit;

                V0_fit(1)=0;
                for k=1:8
                    v_max_fit(k)=alpha_fit(k)*(alpha_fit(k)*tau/(alpha_fit(k)*tauin-V0_fit(k)*(tauin-tau)))^(tau/(tauin-tau));
                    V0_fit(k+1)=V0_fit(k).*exp(-50./tau)+alpha_fit(k)*(tauin./(tauin-tau)).*(exp(-50./tauin)-exp(-50./tau));
                end
                v_max_fit(8)=alpha_fit(8)*(alpha_fit(8)*tau/(alpha_fit(8)*tauin-V0_fit(8)*(tauin-tau)))^(tau/(tauin-tau));
                V0_fit(9)=V0_fit(8).*exp(-550./tau)+alpha_fit(8)*(tauin./(tauin-tau)).*(exp(-550./tauin)-exp(-550./tau));
                v_max_fit(9)=alpha_fit(9)*(alpha_fit(9)*tau/(alpha_fit(9)*tauin-V0_fit(9)*(tauin-tau)))^(tau/(tauin-tau));
                v_difference_fit=v_max_fit-V0_fit;

                b=sum((v_difference_fit-vdifferenceexperiment).^2);
                if b<c
                    c=b;
                    taurec_optimal=taurec_fit;
                    A_optimal=A_fit;
                    U_optimal=U_fit;
                    v_difference_optimal=v_difference_fit;
                end
            end
        end
    end

    
    if  taurec_optimal==taurec_min
        taurec_min=max(taurec_min-num_of_intervals/4*dtaurec,20);
        taurec_max=taurec_max-num_of_intervals/4*dtaurec;
        ee=1;
    end
    if taurec_optimal==taurec_max
        taurec_max=taurec_max+num_of_intervals/2*dtaurec;
        taurec_min=taurec_min+num_of_intervals/2*dtaurec;
        ee1=1;
    end
    if U_optimal==U_min
        U_min=max(0.02,U_min-num_of_intervals/4*dU);
        U_max=U_max-num_of_intervals/4*dU;
        bb=1;
    end
    if U_optimal==U_max
        U_max=U_max+num_of_intervals/4*dU;
        U_min=U_min+num_of_intervals/4*dU;
        bb1=1;
    end
    if A_optimal==A_min
        A_min=max(0.5,A_min-num_of_intervals/4*dA);
        A_max=A_max-num_of_intervals/4*dA;
        cc=1;
    end
    if A_optimal==A_max
        A_max=A_max+num_of_intervals/4*dA;
        A_min=A_min+num_of_intervals/4*dA;
        dd=1;
    end
    if ( ee==0 & ee1==0 & bb==0 & bb1==0 & cc==0 & dd==0 )
        taurec_min=max(20,taurec_optimal-num_of_intervals/4*dtaurec);
        taurec_max=taurec_optimal+num_of_intervals/4*dtaurec;
        dtaurec=(taurec_max-taurec_min)/num_of_intervals;
        U_min=max(0.02,U_optimal-num_of_intervals/4*dU);
        U_max=U_optimal+num_of_intervals/4*dU;
        dU=(U_max-U_min)/num_of_intervals;
        A_min=max(0.5,A_optimal-num_of_intervals/4*dA);
        A_max=A_optimal+num_of_intervals/4*dA;
        dA=(A_max-A_min)/num_of_intervals;
        kk=kk+1
    end
end

% the final values from the fit
cfit=c;   % the overall error, tell us how good the fit is.
Afit=A_optimal;   % the absolute synaptic efficacy.
Ufit=U_optimal;   % the utilization parameter. 
taurecfit=taurec_optimal;   % the recovery time constant. 
taufit=tau_optimal;   % the membrane time constant. 
tauinfit=tauin_optimal;   % the synaptic time constant. 
vdifferencefit=v_difference_optimal;   % the best fitted amplitdues. 

% plotting the experimental and the fitted EPSP amplitudes. 
figure
plot(vdifferenceexperiment,'.-')
hold on
plot(vdifferencefit,'r.-')
box off
axis([0.75 9.25 -0.1 1.95])
xlabel('Input index', 'fontsize', 14)
ylabel('EPSP amplitude [mV]', 'fontsize', 14)
title('EPSP amplitudes: Experiemnt vs. Fit', 'fontsize', 20, 'fontweight', 'b')
legend('Exp.', 'Fit')

save Fit_dynamic_parameters_depression cfit Ufit Afit taurecfit taufit tauinfit vdifferencefit vdifferenceexperiment
