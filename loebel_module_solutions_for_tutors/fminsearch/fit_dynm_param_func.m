function f=fit_dynm_param_func(x,vdifferenceexperiment, tauinfit, taufit)
% This function calculates the mean-squares differnce between the model's
% estimated EPSPs and the experimental EPSPs for a given set of parameters.

taurec_fit=x(1);
U_fit=x(2);
A_fit=x(3);

x_fit(1)=1;
for n=1:7
    x_fit(n+1)=x_fit(n)*(1-U_fit)*exp(-50/taurec_fit)+1-exp(-50/taurec_fit);
end
x_fit(9)=x_fit(8)*(1-U_fit)*exp(-550/taurec_fit)+1-exp(-550/taurec_fit);
alpha_fit=A_fit.*U_fit.*x_fit;

V0_fit(1)=0;
for k=1:7
    v_max_fit(k)=alpha_fit(k)*(alpha_fit(k)*taufit/(alpha_fit(k)*tauinfit-V0_fit(k)*(tauinfit-taufit)))^(taufit/(tauinfit-taufit));
    V0_fit(k+1)=V0_fit(k).*exp(-50./taufit)+alpha_fit(k)*(tauinfit./(tauinfit-taufit)).*(exp(-50./tauinfit)-exp(-50./taufit));
end

v_max_fit(8)=alpha_fit(8)*(alpha_fit(8)*taufit/(alpha_fit(8)*tauinfit-V0_fit(8)*(tauinfit-taufit)))^(taufit/(tauinfit-taufit));
V0_fit(9)=V0_fit(8).*exp(-550./taufit)+alpha_fit(8)*(tauinfit./(tauinfit-taufit)).*(exp(-550./tauinfit)-exp(-550./taufit));
v_max_fit(9)=alpha_fit(9)*(alpha_fit(9)*taufit/(alpha_fit(9)*tauinfit-V0_fit(9)*(tauinfit-taufit)))^(taufit/(tauinfit-taufit));
v_difference_fit=v_max_fit-V0_fit;


f=sum((v_difference_fit-vdifferenceexperiment).^2);






