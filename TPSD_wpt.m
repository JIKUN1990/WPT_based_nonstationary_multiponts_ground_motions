function[dfw,TPSD]=TPSD_wpt(x_input,dt,n)
x=x_input;     
duration=dt*length(x); 
fs=1/dt;
N=length(x_input); 
wpt=wpdec(x_input,n,'dmey');     
%wpt=wpdec(x_input,n,'Sym5');   
% WPT分解
% 实现对节点顺序按照频率递增进行重排序
nodes=get(wpt,'tn');  % Wavelet packet decomposition coefficient
N_cfs=length(nodes);  % Number of wavelet packet coefficients
ord=wpfrqord(nodes);  
nodes_ord=nodes(ord); % Rearranged wavelet coefficients
%%Obtain the wavelet packet coefficients of each frequency band %%%%%%%%%%%%%%%%%%%%%
t=0:2^n*dt:N*dt;
% for i=1:2^n
%     cfs(:,i)=wpcoef(wpt,nodes_ord(i));
% end
for i=1:2^n
    wprc(:,i)=wprcoef(wpt,nodes_ord(i));
end
dfw=(fs/2)/(2^n);
TPSD=wprc.^2/(dfw); % Time-varying PSD
