function[xinter1, xinter2]=FFT_Sa_inter_JK_target_design(x_input,xwpt,dt,n,periods,SAmax)
Lori=length(x_input);
n1=nextpow2(length(x_input));
x1=xwpt;
x_input1=[x_input; zeros(2^n1-length(x_input),1)];
xwpt=[xwpt; zeros(2^n1-length(xwpt),1)];
FAS1=fft(x_input1,2^(n1));
FAS2=fft(x1,2^(n1));
nq=1/dt/2;
df=nq/(2^(n1-1));
% n1=8
% dfw=50/(2^n1);
% Ftar=dfw:dfw:50;
% periods=1./Ftar;
% periods=[0.01 0.013 0.015 0.018 0.0200000000000000,0.0220000000000000,0.0250000000000000,0.0290000000000000,0.0300000000000000,0.0320000000000000,0.0350000000000000,0.0360000000000000,0.0400000000000000,0.0420000000000000,0.0440000000000000,0.0450000000000000,0.0460000000000000,0.0480000000000000,0.0500000000000000,0.0550000000000000,0.0600000000000000,0.0650000000000000,0.0670000000000000,0.0700000000000000,0.0750000000000000,0.0800000000000000,0.0850000000000000,0.0900000000000000,0.0950000000000000,0.100000000000000,0.110000000000000,0.120000000000000,0.130000000000000,0.133000000000000,0.140000000000000,0.150000000000000,0.160000000000000,0.170000000000000,0.180000000000000,0.190000000000000,0.200000000000000,0.220000000000000,0.240000000000000,0.250000000000000,0.260000000000000,0.280000000000000,0.290000000000000,0.300000000000000,0.320000000000000,0.340000000000000,0.350000000000000,0.360000000000000,0.380000000000000,0.400000000000000,0.420000000000000,0.440000000000000,0.450000000000000,0.460000000000000,0.480000000000000,0.500000000000000,0.550000000000000,0.600000000000000,0.650000000000000,0.667000000000000,0.700000000000000,0.750000000000000,0.800000000000000,0.850000000000000,0.900000000000000,0.950000000000000,1,1.10000000000000,1.20000000000000,1.30000000000000,1.40000000000000,1.50000000000000,1.60000000000000,1.70000000000000,1.80000000000000,1.90000000000000,2,2.20000000000000,2.40000000000000,2.50000000000000,2.60000000000000,2.80000000000000,3,3.20000000000000,3.40000000000000,3.50000000000000,3.60000000000000,3.80000000000000,4,4.20000000000000,4.40000000000000,4.60000000000000,4.80000000000000,5,5.50000000000000,6,6.50000000000000,7,7.50000000000000,8,8.50000000000000,9,9.50000000000000,10]
% periods=flip(periods);
T=periods;
% amax=0.2;
% Tg=0.5
mm=length(T);
% for i=1:mm
%     if T(i)<=0.1
%         R(1,i)=5.5*amax*T(i)+0.45*amax;
%     elseif T(i)>0.1 & T(i)<=Tg
%         R(1,i)=amax;
%     elseif T(i)>Tg & T(i)<=5*Tg
%         R(1,i)=(Tg/T(i))^0.9*amax;
%     elseif T(i)>5*Tg & T(i)<=10
%         R(1,i)=[0.2^0.9-0.02*(T(i)-5*Tg)]*amax;
%     end
% end

% SAmax=R*980;

% [SAmax,SVmax,SDmax,location1] = response_spectral_JK_fast(x_input1,[dt:dt:length(x_input1)*dt],periods,0.05);
[SAmax1,SVmax1,SDmax1,location1] = response_spectral_JK_fast(xwpt,[dt:dt:length(xwpt)*dt],periods,0.05);
Ftar1=1./(periods);
Ftar=0:df:nq-df;
% Ftar=df:df:nq;
SAmaxinter=interp1(Ftar1,SAmax,Ftar,'linear');
SAmaxinter2=interp1(Ftar1,SAmax1,Ftar,'linear');
SAmaxinter(isnan(SAmaxinter))=1;
SAmaxinter2(isnan(SAmaxinter2))=1;
scale=SAmaxinter./SAmaxinter2;
scalefft=[scale,fliplr(scale)]

% timetrans=ifft(FAS1);
timetrans1=real(ifft(FAS2.*scalefft'));
xinter1=real(timetrans1(1:length(x_input)));
% PGAscale=max(abs(x_input))./max(abs(xinter1));
% Egscale=max(cumsum(x_input.^2))/max(cumsum(xinter1.^2));
% xinter1=sqrt(Egscale).*xinter1;

[dfw,cfs,nodes_ord]=CFS_wpt(xwpt,dt,n);
[mm,nn]=size(cfs);
kk=floor(length(x_input)/mm)
Ia=cumsum(x_input(1:kk*mm).^2)
Ia1=cumsum(xinter1(1:kk*mm).^2)

for i=1:mm
 scale2(i)= sqrt(Ia(kk+(i-1)*kk)-Ia(1+(i-1)*kk))./sqrt(Ia1(kk+(i-1)*kk)-Ia1(1+(i-1)*kk));
end
for i=1:mm
xinter2(1+(i-1)*kk:kk+(i-1)*kk)=xinter1(1+(i-1)*kk:kk+(i-1)*kk).*scale2(i);
end
xinter2=[xinter2';zeros(length(x_input)-length(xinter2),1)];
