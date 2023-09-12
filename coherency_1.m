function [f1H,r12]=coherency_1(x_input,x1,x2,dt)
dAH_1=x1;
dAH_2=x2;


Fs=1/dt;            
num1H=length(dAH_1);              
N1H=2^nextpow2(num1H);          
dAH_1=[dAH_1;zeros(N1H-length(dAH_1),1)];
ff1H=(0:1:N1H-1)/N1H*Fs; 
f1H=ff1H';
y1H=fft(dAH_1,N1H)*dt;      
fH=(0:1:N1H-1)/N1H*Fs;    
S1=conj(y1H).*y1H;        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


num2H=length(dAH_2);     
N2H=2^nextpow2(num2H);           
dAH_2=[dAH_2;zeros(N1H-length(dAH_2),1)];
f2H=(0:1:N2H-1)/N2H*Fs; 
y2H=fft(dAH_2,N2H)*dt;     
fH=(0:1:N2H-1)/N2H*Fs;  


S1=conj(y1H).*y1H;
S2=conj(y2H).*y2H;
S12=conj(y1H).*y2H;

 L=50;
 w = hamming(L);
 w1=w/(sum(w));
 
 S1_G = conv2(S1,w1,'same');
 S2_G = conv2(S2,w1,'same');
 S12_G = conv2(S12,w1,'same');
 
yyy=abs(S12_G)./sqrt(S1_G.*S2_G);
r12=yyy;
% n=8
% df1=(Fs/2)./(2^n)
% % Ttar=df1:df1:100;
% Ttar=df1:df1:(1/dt);
% % cfs_w=dfw-0.5*dfw:dfw:(2^n)*dfw-0.5*dfw
% Nq=1/(2*dt)
% Ttar=df1:df1:Nq*2
% r12 = interp1(f1H,yyy,Ttar,'linear');
% r12(end)=r12(end-1);
% plot(Ttar,yyy2)




% subplot(4,1,1)
% plot(dt:dt:dt*length(dAH_1),dAH_1)
% hold on
% plot(dt:dt:dt*length(x_input),x_input,'g','lineWidth',3)
% hold on
% subplot(4,1,2)
% plot(dt:dt:dt*length(dAH_2),dAH_2)
% hold on
% plot(dt:dt:dt*length(dAH_1),dAH_1)
% subplot(4,1,3)
% plot(f1H,yyy);
% xlim([0 10]);

% subplot(4,1,4)
% plot(fH,S1);
% hold on
% plot(fH,S2);
% xlim([0 10]);
end
