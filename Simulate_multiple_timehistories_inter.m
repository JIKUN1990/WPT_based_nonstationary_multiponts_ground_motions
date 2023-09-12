
function[y12,y13,y23,x1,x2,x3]=Simulate_multiple_timehistories_inter(TPSD,N,dt,dfw,R,d12,d23,d13,v_app,indexcase)
% n=length(TPSD(1,:));
rand1=R.rand1;
rand2=R.rand2;
rand3=R.rand3;

% cfs_w=dfw:dfw:2^n*dfw;
% w=cfs_w;
% w1=0.1:0.1:N*0.1;
dw=0.05
w1=dw:dw:N*dw;
cfs_w=w1./(2*pi);
w=cfs_w;
thetaw=5210*(1+(w/1.09).^2.78).^(-0.5);
y12=0.736*exp(-2*d12./(0.147*thetaw).*(1-0.736+0.147*0.736))+(1-0.736).*exp(-2*d12./(thetaw).*(1-0.736+0.147*0.736));
y13=0.736*exp(-2*d13./(0.147*thetaw).*(1-0.736+0.147*0.736))+(1-0.736).*exp(-2*d13./(thetaw).*(1-0.736+0.147*0.736));
y23=0.736*exp(-2*d23./(0.147*thetaw).*(1-0.736+0.147*0.736))+(1-0.736).*exp(-2*d23./(thetaw).*(1-0.736+0.147*0.736));

y21=y12;y31=y13;y32=y23;

for i=1:N
    B=  [1 y12(i) y13(i);y21(i) 1 y23(i); y31(i) y32(i) 1];
    beta(:,:,i)=chol(B)';
end
% v_app=1000
phase11=zeros([1,length(cfs_w)]);
phase12=-2*pi*cfs_w.*d12./v_app;
phase13=-2*pi*cfs_w.*d13./v_app;
phase21=-phase12;
phase22=zeros([1,length(cfs_w)]);
phase23=-2*pi*cfs_w.*d23./v_app;
phase32=-phase23;phase31=-phase13;
phase33=zeros([1,length(cfs_w)]);

beta11=reshape(beta(1,1,:),[1,length(cfs_w)]);
beta12=reshape(beta(1,2,:),[1,length(cfs_w)]);
beta13=reshape(beta(1,3,:),[1,length(cfs_w)]);
beta21=reshape(beta(2,1,:),[1,length(cfs_w)]);
beta22=reshape(beta(2,2,:),[1,length(cfs_w)]);
% beta22=ones([1,length(cfs_w)]);
beta23=reshape(beta(2,3,:),[1,length(cfs_w)]);
beta31=reshape(beta(3,1,:),[1,length(cfs_w)]);
beta32=reshape(beta(3,2,:),[1,length(cfs_w)]);
beta33=reshape(beta(3,3,:),[1,length(cfs_w)]);
% beta33=ones([1,length(cfs_w)]);

% %%
% for i=1:length(TPSD)
%     i;
%     step=i;
%     s=(TPSD(i,:));
%     S11=s; S22=s;  S33=s;

% s=TPSD;
S11=TPSD.TPSD1; 
S22=TPSD.TPSD2;  
S33=TPSD.TPSD3;

T=dt:dt:dt*length(S11(:,1));
T=T';
if indexcase==1
Time11=2*sqrt(0.1).*sqrt(S11).*beta11.*cos(2*pi*cfs_w.*T-phase11+rand1);
Time21=2*sqrt(0.1).*sqrt(S22).*beta21.*cos(2*pi*cfs_w.*T-phase21+rand1);
Time22=2*sqrt(0.1).*sqrt(S22).*beta22.*cos(2*pi*cfs_w.*T-phase22+rand2);
Time31=2*sqrt(0.1).*sqrt(S33).*beta31.*cos(2*pi*cfs_w.*T-phase31+rand1);
Time32=2*sqrt(0.1).*sqrt(S33).*beta32.*cos(2*pi*cfs_w.*T-phase32+rand2);
Time33=2*sqrt(0.1).*sqrt(S33).*beta33.*cos(2*pi*cfs_w.*T-phase33+rand3);
elseif indexcase==0
Time11=2*sqrt(0.1).*sqrt(S11).*beta11.*cos(2*pi*cfs_w.*T+rand1);
Time21=2*sqrt(0.1).*sqrt(S22).*beta21.*cos(2*pi*cfs_w.*T+rand1);
Time22=2*sqrt(0.1).*sqrt(S22).*beta22.*cos(2*pi*cfs_w.*T+rand2);
Time31=2*sqrt(0.1).*sqrt(S33).*beta31.*cos(2*pi*cfs_w.*T+rand1);
Time32=2*sqrt(0.1).*sqrt(S33).*beta32.*cos(2*pi*cfs_w.*T+rand2);
Time33=2*sqrt(0.1).*sqrt(S33).*beta33.*cos(2*pi*cfs_w.*T+rand3);
end

% Time11=2*sqrt(dw).*sqrt(S11).*beta11.*cos(2*pi*cfs_w.*T+rand1);
% Time21=2*sqrt(dw).*sqrt(S22).*beta21.*cos(2*pi*cfs_w.*T+rand1);
% Time22=2*sqrt(dw).*sqrt(S22).*beta22.*cos(2*pi*cfs_w.*T+rand2);
% Time31=2*sqrt(dw).*sqrt(S33).*beta31.*cos(2*pi*cfs_w.*T+rand1);
% Time32=2*sqrt(dw).*sqrt(S33).*beta32.*cos(2*pi*cfs_w.*T+rand2);
% Time33=2*sqrt(dw).*sqrt(S33).*beta33.*cos(2*pi*cfs_w.*T+rand3);
Point3=Time31+Time32+Time33;
Point2=Time21+Time22;
Point1=Time11;
x1=sum(Point1')'; % 将不同频率层级的信号组合起来得到最终记录
x2=sum(Point2')'; %
x3=sum(Point3')'; %