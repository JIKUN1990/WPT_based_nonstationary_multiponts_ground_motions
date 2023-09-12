function[dfw,cfs,nodes_ord]=CFS_wpt(x_input,dt,n)
x=x_input;        %   查看频谱范围
% duration=dt*length(x); % 记录持时
fs=1/dt;
N=length(x_input); %采样点个数
% n=8;
wpt=wpdec(x_input,n,'dmey');      %用meyr小波进行7层小波包分解%%%%%%%%%%%%%%%%%%%%%%
% WPT分解
% 实现对节点顺序按照频率递增进行重排序
nodes=get(wpt,'tn');  %小波包分解系数
N_cfs=length(nodes);  %小波包系 数个数
ord=wpfrqord(nodes);  %小波包系数重排，ord是重排后小波包系数索引构成的矩阵　如3层分解的[1;2;4;3;7;8;6;5]
nodes_ord=nodes(ord); %重排后的小波系数
%%求取小波包分解的各个频段的小波包系数%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=0:2^n*dt:N*dt;
for i=1:2^n
    cfs(:,i)=wpcoef(wpt,nodes_ord(i));
end
%%时频功率谱
% for i=1:2^n
%     wprc(:,i)=wprcoef(wpt,nodes_ord(i));
% end
dfw=(fs/2)/(2^n);
% TPSD=wprc.^2/(dfw); % 时频功率谱
