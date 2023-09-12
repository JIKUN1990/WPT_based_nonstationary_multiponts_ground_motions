function [kk]=plotFig_results(x1,x2,x3,dt,n,y12,y13,y23,T_target,SAtarget)
% FAS
x_input=x1;
dfw=((1/dt)./2)/(2^n);
[mm,~]=size(x1);
if mm==1
    x1=x1';
    x2=x2';
    x3=x3';
end
dfw=(1/dt/2)/(2^n);
subplot(4,2,7);
Ftar=dfw:dfw:(1/dt/2);
periods=1./Ftar;
% TT=periods
TT=[0.01 0.012 0.013 0.015 0.018 0.019 0.0200000000000000,0.0220000000000000,0.0250000000000000,0.0290000000000000,0.0300000000000000,0.0320000000000000,0.0350000000000000,0.0360000000000000,0.0400000000000000,0.0420000000000000,0.0440000000000000,0.0450000000000000,0.0460000000000000,0.0480000000000000,0.0500000000000000,0.0550000000000000,0.0600000000000000,0.0650000000000000,0.0670000000000000,0.0700000000000000,0.0750000000000000,0.0800000000000000,0.0850000000000000,0.0900000000000000,0.0950000000000000,0.100000000000000,0.110000000000000,0.120000000000000,0.130000000000000,0.133000000000000,0.140000000000000,0.150000000000000,0.160000000000000,0.170000000000000,0.180000000000000,0.190000000000000,0.200000000000000,0.220000000000000,0.240000000000000,0.250000000000000,0.260000000000000,0.280000000000000,0.290000000000000,0.300000000000000,0.320000000000000,0.340000000000000,0.350000000000000,0.360000000000000,0.380000000000000,0.400000000000000,0.420000000000000,0.440000000000000,0.450000000000000,0.460000000000000,0.480000000000000,0.500000000000000,0.550000000000000,0.600000000000000,0.650000000000000,0.667000000000000,0.700000000000000,0.750000000000000,0.800000000000000,0.850000000000000,0.900000000000000,0.950000000000000,1,1.10000000000000,1.20000000000000,1.30000000000000,1.40000000000000,1.50000000000000,1.60000000000000,1.70000000000000,1.80000000000000,1.90000000000000,2,2.20000000000000,2.40000000000000,2.50000000000000,2.60000000000000,2.80000000000000,3,3.20000000000000,3.40000000000000,3.50000000000000,3.60000000000000,3.80000000000000,4,4.20000000000000,4.40000000000000,4.60000000000000,4.80000000000000,5,5.50000000000000,6,6.50000000000000,7,7.50000000000000,8,8.50000000000000,9,9.50000000000000,10];
% [SAmax,SVmax,SDmax ] = response_spectral_JK_fast(x_input,[0.01:0.01:length(x_input)*0.01],TT,0.05)
% plot(TT,SAmax,'lineWidth',1)
% hold on
[SAmax1,SVmax,SDmax ] = response_spectral_JK_fast(x1,[0.01:0.01:length(x_input)*0.01],TT,0.05);
plot(TT,SAmax1,'lineWidth',1);
hold on
[SAmax2,SVmax,SDmax ] = response_spectral_JK_fast(x2,[0.01:0.01:length(x_input)*0.01],TT,0.05);
plot(TT,SAmax2,'lineWidth',1);
hold on
[SAmax3,SVmax,SDmax ] = response_spectral_JK_fast(x3,[0.01:0.01:length(x_input)*0.01],TT,0.05);
plot(TT,SAmax3,'lineWidth',1);
hold on
plot(T_target,SAtarget,'lineWidth',3);
legend('point1','point2','point3','target spectrum','location','northeast');

% figure(2)
% subplot(3,3,5)
% plot(cumsum(x1.^2))
% hold on
% plot(cumsum(x2.^2))
% hold on
% plot(cumsum(x3.^2))
% hold on
% plot(cumsum(x_input.^2))
subplot(4,2,2)
[f12,r12]=coherency_1(x_input,x1,x2,dt);
plot(f12,r12);
xlim([0 10])
hold on
plot(dfw:dfw:2^n.*dfw,y12,'lineWidth',3)
title('coherency between point1 and point2')

subplot(4,2,4)
[f13,r13]=coherency_1(x_input,x1,x3,dt);
plot(f13,r13);
xlim([0 10])
hold on
plot(dfw:dfw:2^n.*dfw,y13,'lineWidth',3)
title('coherency between point1 and point3')
subplot(4,2,6)
[f23,r23]=coherency_1(x_input,x2,x3,dt);
plot(f23,r23);
xlim([0 10])
hold on
plot(dfw:dfw:2^n.*dfw,y23,'lineWidth',3)
title('coherency between point2 and point3')
xlabel('Frequency(Hz)')


subplot(4,2,1)
% plot(dt:dt:length(x_input)*dt,x_input,'r','lineWidth',1)
% hold on
plot(dt:dt:length(x1)*dt,x1,'k','lineWidth',0.5)
hold on
ylabel('gal')
legend('point1')
subplot(4,2,3)
% plot(dt:dt:length(x_input)*dt,x_input,'r','lineWidth',1)
% hold on
plot(dt:dt:length(x1)*dt,x1,'k','lineWidth',0.5)
ylabel('gal')
legend('point2')
subplot(4,2,5)
% plot(dt:dt:length(x_input)*dt,x_input,'r','lineWidth',1)
% hold on
plot(dt:dt:length(x1)*dt,x1,'k','lineWidth',0.5)
xlabel('Time(s)')
ylabel('gal')
legend('point3')
set(gcf,'position',[100 100 1000 1000]);











% subplot(3,3,3)
% [f12,r12]=coherency_2(x_input,x1,x2,dt);
% plot(f12,r12);
% hold on
% plot(dfw:dfw:2^n.*dfw,y12,'lineWidth',3)
% xlim([0 10])
% 
% subplot(3,3,6)
% [f13,r13]=coherency_2(x_input,x1,x3,dt);
% plot(f13,r13);
% xlim([0 10])
% hold on
% plot(dfw:dfw:2^n.*dfw,y13,'lineWidth',3)
% 
% subplot(3,3,9)
% [f23,r23]=coherency_2(x_input,x2,x3,dt);
% plot(f23,r23);
% xlim([0 10])
% hold on
% plot(dfw:dfw:2^n.*dfw,y23,'lineWidth',3)
