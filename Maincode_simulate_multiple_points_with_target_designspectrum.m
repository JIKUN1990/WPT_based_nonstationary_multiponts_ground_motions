% Program for generating multiplle spatially correlated non-stationary 
% ground motions by wavelet-packet algorithm and spectral representation method 
% coded by Associate Professor Kun Ji(Hohai University) 2023/09/06
% jikun@hhu.edu.cn
% Please cite as
% " Ji, K., Cao, X., Wang, S., & Wang, S. (2023). 
% Efficient simulation of spatially correlated non-stationary ground motions by wavelet-packet algorithm and spectral representation method 
% Earthquake Engineering and Engineering Vibration (under review)
clc
clear
%% Input parameters
Mw=7.25; % Magnitude parameter; (earthquake scenario disaggregation);
Rrup=18.7; % Rupture distance parameter; (earthquake scenario disaggregation);
Vs30=760; %Soil condition of target site;(m/s)
ncase=3; %number of generated seed recordings；
% Input parameters for generating ground motions at multiple points 
n=8; % WPT decomposition level, setting default values as 8;
num_case=10; % Number of realizations for generated spatially varying ground motions;
d12=200; %unit:m
d23=200; %unit:m
d13=400; %unit:m
Dis_points=[d12;d23;d13]; % Distance between target sites;
v_app=2000 % apparent wave velocity; unit:m/s
periods=flip([0.010 0.013 0.015 0.018 0.020 0.022 0.025 0.029 0.030 0.032 0.035 0.036 0.040 0.042 0.044 0.045 0.046 	0.048 	0.050 	0.055 	0.060 	0.065 	0.067 	0.070 	0.075 	0.080 	0.085 	0.090 	0.095 	0.100 	0.110 	0.120 	0.130 	0.133 	0.140 	0.150 	0.160 	0.170 	0.180 	0.190 	0.200 	0.220 	0.240 	0.250 	0.260 	0.280 	0.290 	0.300 	0.320 	0.340 	0.350 	0.360 	0.380 	0.400 	0.420 	0.440 	0.450 	0.460 	0.480 	0.500 	0.550 	0.600 	0.650 	0.667 	0.700 	0.750 	0.800 	0.850 	0.900 	0.950 	1.000 	1.100 	1.200 	1.300 	1.400 	1.500 	1.600 	1.700 	1.800 	1.900 	2.000 	2.200 	2.400 	2.500 	2.600 	2.800 	3.000 	3.200 	3.400 	3.500 	3.600 	3.800 	4.000 	4.200 	4.400 	4.600 	4.800 	5.000 	5.500 	6.000 	6.500 	7.000 	7.500 	8.000 	8.500 	9.000 	9.500 	10.000]);
%% Load or define the target spectrum
load DesignSpectrum.mat % T refers to the periods; SAtarget refers to the target design spectrum, with unit g.
default_duration=20; 
P=0.85;
w=(1./periods);
temp=-pi*log(P)./(default_duration*w);
PSD=-0.05./(pi*w).*((SAtarget./980).^2)./log(temp);
E_target=trapz(w,PSD*(default_duration+10)); 
totalEnergytarget=log(E_target); % Calculate the total target energy for design response spectrum
%% genertaion of seed ground motions with total energy as conditional target value
% Program for generating artificial ground motions by stochastic ground
% motion model using wavelet packets
% coded by Yoshi Yamamoto (Stanford University) 10/31/2010
% yama4423@stanford.edu
% The program is modified on the basis of Yamamoto's code, introducing
% total engergy as conditioanl target value to derive other 13 parameters
% of wavelet packets coefficients.
Seedgroundmotion=main_sgm_conditional(Mw, Rrup, Vs30, totalEnergytarget,ncase);
seed_gr_index=1; % selected seed ground motions in n simulated results
T_input=Seedgroundmotion{seed_gr_index,1}(:,1);
x_input=Seedgroundmotion{seed_gr_index,1}(:,2).*980; 
dt=T_input(2)-T_input(1); % time-interal of seed ground motions; 
%% WPT-based simulation program generating spatially varying ground motions
wcutoff=25; % Cutoff frequency
waveindex=0;
% waveindex=1; considering the phase change due to wave propagation effect 
% waveindex=0; without considering the phase change due to wave propagation effect 
[y12,y13,y23,Result,RC12,RC23,RC13]=multipoint_simulate(Dis_points,v_app,x_input,wcutoff,dt,n,num_case,waveindex)
% y12,y13,y23; target empirical coherence model
% cfs_w=dfw:dfw:2^n*dfw;
% w1=cfs_w;
% thetaw=5210*(1+(w1/1.09).^2.78).^(-0.5);
% y12=0.736*exp(-2*d12./(0.147*thetaw).*(1-0.736+0.147*0.736))+(1-0.736).*exp(-2*d12./(thetaw).*(1-0.736+0.147*0.736));
% y13=0.736*exp(-2*d13./(0.147*thetaw).*(1-0.736+0.147*0.736))+(1-0.736).*exp(-2*d13./(thetaw).*(1-0.736+0.147*0.736));
% y23=0.736*exp(-2*d23./(0.147*thetaw).*(1-0.736+0.147*0.736))+(1-0.736).*exp(-2*d23./(thetaw).*(1-0.736+0.147*0.736));
% x11,x22,x33：coherency between generated ground motions at multi-points  
% x1,x2,x3; generated ground motions at multi-points;  
% RC12,RC23,RC13; coherency between generated ground motions at multi-points  
%% Compare the coherency loss with target coherence model 
dfw=(1/dt)/(2*2^n);
x1=Result.Timehistories{1,1}(:,1);
x2=Result.Timehistories{1,1}(:,2);
x3=Result.Timehistories{1,1}(:,3);
[f12,~]=coherency_1(x_input,x2,x3,dt);
[f13,~]=coherency_1(x_input,x2,x3,dt);
[f23,~]=coherency_1(x_input,x2,x3,dt);
cfs_w=dfw:dfw:2^n*dfw;
subplot(4,1,1)
plot(f12,RC12,'g')
xlim([0 10]);
hold on
plot(f12,mean(RC12'),'lineWidth',3);
hold on
plot(cfs_w,y12,'lineWidth',3)
title('coherency between point1 and point2')

subplot(3,1,2)
plot(f13,RC13,'g')
xlim([0 10]);
hold on
plot(f13,mean(RC13'),'lineWidth',3);
hold on
plot(cfs_w,y13,'lineWidth',3)
title('coherency between point1 and point3')

subplot(3,1,3)
plot(f23,RC23,'g')
xlim([0 10]);
hold on
plot(f23,mean(RC23'),'lineWidth',3);
hold on
plot(cfs_w,y23,'lineWidth',3)
title('coherency between point2 and point3')
%% Spectrum matching
pickindex=1 %pick statisfactory simulated results in multiple realizations
x1=Result.Timehistories{pickindex,1}(:,1);
x2=Result.Timehistories{pickindex,1}(:,2);
x3=Result.Timehistories{pickindex,1}(:,3);
for i=1:5
    i;
    if i==1
    [xinter1 xinter2]=FFT_Sa_inter_JK_target_design(x_input,x1,dt,n,periods,SAtarget);
    else 
     [xinter1 xinter2]=FFT_Sa_inter_JK_target_design(x_input,xinter1,dt,n,periods,SAtarget);
    end
end
x11=xinter1;
for i=1:5
    i;
    if i==1
    [xinter1 xinter2]=FFT_Sa_inter_JK_target_design(x_input,x2,dt,n,periods,SAtarget);
    else 
     [xinter1 xinter2]=FFT_Sa_inter_JK_target_design(x_input,xinter1,dt,n,periods,SAtarget);
    end
end
x22=xinter1;
for i=1:5
    i;
    if i==1
    [xinter1 xinter2]=FFT_Sa_inter_JK_target_design(x_input,x3,dt,n,periods,SAtarget);
    else 
     [xinter1 xinter2]=FFT_Sa_inter_JK_target_design(x_input,xinter1,dt,n,periods,SAtarget);
    end
end
x33=xinter1;
%% Plot results
plotFig_results(x11,x22,x33,dt,n,y12,y13,y23,periods,SAtarget)
%%
saveas(gcf,'Generated multiple ground motions result.png')


