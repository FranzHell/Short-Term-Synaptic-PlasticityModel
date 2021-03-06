spiketimes = [100:50:450 1000]
tau_mem =32;
tau_in = 1.8;
A = 144;
u = 0.26;
tau_rec =  1000;
dt = 0.25;

x(1) = 1;
y(1) = 0
V(1) = 0;

t = 0;
tmax = 1200;

ind = 1;
kk = 1;
while t< tmax
    %if spike arrives
    if t>= spiketimes(kk) & t<spiketimes(kk)+dt
        x(ind+1) =  x(ind) + dt*((1-x(ind))/tau_rec)-(u*x(ind)/dt);
        y(ind+1) =  y(ind) + dt*(-y(ind)/tau_in)+(u*x(ind)/dt);
        if kk <9
            kk =  kk+1;
        end
        else
      
        x(ind+1) =  x(ind) + dt*((1-x(ind))/tau_rec);
        y(ind+1) =  y(ind) + dt*(-y(ind)/tau_in);
        end
        
        V(ind+1) = V(ind) + dt/tau_mem*(-V(ind)+A*y(ind));
        
        ind = ind+1;
        t = t+dt;
end

%     figure
%     plot(t:dt:tmax,x,'k.-')
%     title('resources')
%     
%     figure
%     plot(t:dt:tmax,y,'k.-')
%     
%     figure
%     plot(t:dt:tmax,V)
%     
%     load Data_depressing
%     plot(mean(V_depressing))


figure(1)
plot(0:0.25:tmax,x,'k.-')
xlabel('Time [ms]', 'Fontsize', 16)
ylabel('x', 'Fontsize', 16)
title('Fraction of available resources', 'Fontsize', 20, 'Fontweight', 'B')
box off
axis ([-0.5 1220 -0.1 1.1])

figure(2)
plot(0:0.25:tmax,y,'k.-')
xlabel('Time [ms]', 'Fontsize', 16)
ylabel('y', 'Fontsize', 16)
title('Fraction of active resources', 'Fontsize', 20, 'Fontweight', 'B')
box off
axis ([-0.5 1220 -0.1 1.1])

figure(3)
plot(0:0.25:tmax,V,'k.-')
xlabel('Time [ms]', 'Fontsize', 16)
ylabel('V [mV]', 'Fontsize', 16)
title('Simulated membrane potential', 'Fontsize', 20, 'Fontweight', 'B')
box off
axis ([-0.5 1220 -0.1 2.1])

load Data_depressing
figure(4)
plot(0.25:0.25:1200, mean(V_depressing), 'b.-')
xlabel('Time [ms]', 'Fontsize', 16)
ylabel('V [mV]', 'Fontsize', 16)
title('Measured membrane potential', 'Fontsize', 20, 'Fontweight', 'B')
box off
axis ([-0.5 1220 -0.1 2.1])




