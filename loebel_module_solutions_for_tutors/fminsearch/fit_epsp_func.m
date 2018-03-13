function f=fit_epsp_func(x,y,t,dt)

tauin=x(1);
tau=x(2);
A=x(3);

v=A.*(tauin/(tauin-tau)).*(exp(-(t.*dt)./tauin)-exp(-(t.*dt)./tau));
f=sum((y-v).^2);