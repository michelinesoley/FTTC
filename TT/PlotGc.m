function [memo,xi,xic,t] = PlotGc(G0_tt,G_tt,Gc_tt,xi,xic,d,k,t,tau,dx,nsteps)
figure(7)
% Survival amplitude
xi(k)=0;
t(k)=(k-1)*tau;   % current time
xi(k)=dot(conj(G0_tt),G_tt)*dx^d;
xic(k)=0;
xic(k)=dot(conj(G0_tt),Gc_tt)*dx^d;

subplot(1,2,1);
hold on;
plot(t(1:k)*2.418E-2,real(xi(1:k)),'k:','MarkerSize',15,'linewidth',2);
scatter(t(1:k)*2.418E-2,real(xic(1:k)),250,'filled','b');
title('Real Part','fontsize',18);
ylabel('Real Part','fontsize',18);
xlabel('Time [fs]','fontsize',18);
xlim([0 nsteps*tau*2.418E-2]);

subplot(1,2,2);
hold on;
plot(t(1:k)*2.418E-2,imag(xi(1:k)),'k:','MarkerSize',15,'linewidth',2);
scatter(t(1:k)*2.418E-2,imag(xic(1:k)),250,'filled','b');
title('Imaginary Part','fontsize',18);
ylabel('Imaginary Part, rad','fontsize',18);
xlabel('Time [fs]','fontsize',18);
legend('SOFT','Chebyshev')
title('Imaginary Part','fontsize',18);
ylabel('Imaginary Part','fontsize',18);
xlabel('Time, fs','fontsize',18);
xlim([0 nsteps*tau*2.418E-2]);
hold off;
memo=mem(G_tt); % memory

end

