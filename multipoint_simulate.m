function[y12,y13,y23,Result,RC12,RC23,RC13]= multipoint_simulate(Dis_points,v_app,x_input,wcutoff,dt,n,num_case,index)

dw=0.05; % the frequency step after linear interpolation
TPSD=[];
[dfw TPSD0]=TPSD_wpt(x_input,dt,n); % Generating Time-varying PSD
clf
cfs_w=dfw:dfw:2^n*dfw;
W=2*pi*cfs_w;
w=dw:dw:W(end);
TPSD1=zeros(length(x_input),length(w));
for i=1:length(TPSD0(:,1))
    i;
    TPSD1(i,:)=interp1(W,TPSD0(i,:),w,'linear'); %linear interpolation in frequency domain
end
TPSD1(isnan(TPSD1))=0;
% wcutoff=25
N1=floor(wcutoff*(2*pi)/dw);
N2=floor(wcutoff/dfw);
TPSD1(:,N1:end)=0;
TPSD0(:,N2:end)=0;
A=sum(sum(TPSD0))%Total Engergy before interpolation
B=sum(sum(TPSD1))%Total Engergy after interpolation
TPSD1=TPSD1*(A/B);

close all
% load TPSD1
for Ntimes=1:num_case 
Ntimes;
N=length(TPSD1(1,:));
d12=Dis_points(1);d23=Dis_points(2);d13=Dis_points(3);
% v_app=2000;
TPSD.TPSD1=TPSD1;
TPSD.TPSD2=TPSD1;
TPSD.TPSD3=TPSD1;

interation=1
for NN=1:interation
    
for i=1:N
    R.rand1(i)=unifrnd(0,2*pi,1);
    R.rand2(i)=unifrnd(0,2*pi,1);
    R.rand3(i)=unifrnd(0,2*pi,1);
end
[y12,y13,y23,x1,x2,x3]=Simulate_multiple_timehistories_inter(TPSD,N,dt,dfw,R,d12,d23,d13,v_app,index);
k1=max(cumsum(x_input.^2))/max(cumsum(x1.^2));
k2=max(cumsum(x_input.^2))/max(cumsum(x2.^2));
k3=max(cumsum(x_input.^2))/max(cumsum(x3.^2));
% x1=sqrt(k1).*x1;
% x2=sqrt(k2).*x2;
% x3=sqrt(k3).*x3;

% k1=max(abs(x_input))/max(abs(x1));
% k2=max(abs(x_input))/max(abs(x2));
% k3=max(abs(x_input))/max(abs(x3));


%[S1_G,S2_G,S3_G,TPSD]=Sa_ScaleFunction_2(TPSD1,N,x1,x2,x3,x_input);
end

K_test(Ntimes,1)=sqrt(k1);K_test(Ntimes,2)=sqrt(k2);K_test(Ntimes,3)=sqrt(k3);
TT=[0.0100000000000000,0.0200000000000000,0.0220000000000000,0.0250000000000000,0.0290000000000000,0.0300000000000000,0.0320000000000000,0.0350000000000000,0.0360000000000000,0.0400000000000000,0.0420000000000000,0.0440000000000000,0.0450000000000000,0.0460000000000000,0.0480000000000000,0.0500000000000000,0.0550000000000000,0.0600000000000000,0.0650000000000000,0.0670000000000000,0.0700000000000000,0.0750000000000000,0.0800000000000000,0.0850000000000000,0.0900000000000000,0.0950000000000000,0.100000000000000,0.110000000000000,0.120000000000000,0.130000000000000,0.133000000000000,0.140000000000000,0.150000000000000,0.160000000000000,0.170000000000000,0.180000000000000,0.190000000000000,0.200000000000000,0.220000000000000,0.240000000000000,0.250000000000000,0.260000000000000,0.280000000000000,0.290000000000000,0.300000000000000,0.320000000000000,0.340000000000000,0.350000000000000,0.360000000000000,0.380000000000000,0.400000000000000,0.420000000000000,0.440000000000000,0.450000000000000,0.460000000000000,0.480000000000000,0.500000000000000,0.550000000000000,0.600000000000000,0.650000000000000,0.667000000000000,0.700000000000000,0.750000000000000,0.800000000000000,0.850000000000000,0.900000000000000,0.950000000000000,1,1.10000000000000,1.20000000000000,1.30000000000000,1.40000000000000,1.50000000000000,1.60000000000000,1.70000000000000,1.80000000000000,1.90000000000000,2,2.20000000000000,2.40000000000000,2.50000000000000,2.60000000000000,2.80000000000000,3,3.20000000000000,3.40000000000000,3.50000000000000,3.60000000000000,3.80000000000000,4,4.20000000000000,4.40000000000000,4.60000000000000,4.80000000000000,5,5.50000000000000,6,6.50000000000000,7,7.50000000000000,8,8.50000000000000,9,9.50000000000000,10];


% figure(1)
[f12,r12]=coherency_1(x_input,x1,x2,dt);
% hold on
% plot(cfs_w,y12)
% figure(2)
[f13,r13]=coherency_1(x_input,x1,x3,dt);
% hold on
% plot(cfs_w,y13)
% figure(3)
[f23,r23]=coherency_1(x_input,x2,x3,dt);
% hold on
% plot(cfs_w,y23)
Result.Timehistories{Ntimes,1}=[x1 x2 x3];
Result.Coherency{Ntimes,1}=[r12,r13,r23];
RC12(:,Ntimes)=r12;
RC23(:,Ntimes)=r23;
RC13(:,Ntimes)=r13;
end

% close all
cfs_w=dfw:dfw:2^n*dfw;
w1=cfs_w;
thetaw=5210*(1+(w1/1.09).^2.78).^(-0.5);
y12=0.736*exp(-2*d12./(0.147*thetaw).*(1-0.736+0.147*0.736))+(1-0.736).*exp(-2*d12./(thetaw).*(1-0.736+0.147*0.736));
y13=0.736*exp(-2*d13./(0.147*thetaw).*(1-0.736+0.147*0.736))+(1-0.736).*exp(-2*d13./(thetaw).*(1-0.736+0.147*0.736));
y23=0.736*exp(-2*d23./(0.147*thetaw).*(1-0.736+0.147*0.736))+(1-0.736).*exp(-2*d23./(thetaw).*(1-0.736+0.147*0.736));


% for ii=1:Ntimes
% x11=Result.Timehistories{ii,1}(:,1);
% x22=Result.Timehistories{ii,1}(:,2);
% x33=Result.Timehistories{ii,1}(:,3);
% [t01 t02 Ds0]=Ds595_JK(x_input,dt)
% [t11 t12 Ds1]=Ds595_JK(x11,dt)
% [t21 t22 Ds2]=Ds595_JK(x22,dt)
% [t31 t32 Ds3]=Ds595_JK(x33,dt)
% PGA1=max(abs(x_input))
% PGA2=max(abs(x11))
% PGA3=max(abs(x22))
% PGA4=max(abs(x33))
% Ds(ii,:)=[Ds0 Ds1 Ds2 Ds3]
% PGA(ii,:)=[PGA1 PGA2 PGA3 PGA4]
end

