
clc;
clear; 
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Parameter settings %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PAR.FN=4;   % variable nodes (VN), number of data layers
PAR.VN=6;   % function nodes (FN), number of physical resources
PAR.d_f=3;      % Each FN is connected to 3 VNs
PAR.d_v=2;      % Each VN is connected to 2 FNs 
PAR.M=4;        % Number of codeword in each codebook
PAR.Data_length=200;    % Number of total data
PAR.N_iter= 5;  % Number of iterations in decoding
%PAR.EbNo = 10;
PAR.EbNo = [1:20];
PAR.CB=zeros(PAR.M,PAR.FN,PAR.VN);      %Codebooks
%{
PAR.CB(:,:,1)=[  0,0,0,0;
            -0.1815-0.1318i,-0.6351-0.4615i,0.6351+0.4615i,0.1815+0.1318i;
            0,0,0,0;
            0.7851,-0.2243,0.2243,-0.7851];
 
PAR.CB(:,:,2)=[  0.7851,-0.2243,0.2243,-0.7851;
            0,0,0,0;
            -0.1815-0.1318i,-0.6351-0.4615i,0.6351+0.4615i,0.1815+0.1318i;
            0,0,0,0];
         
PAR.CB(:,:,3)=[  -0.6351+0.4615i,0.1815-0.1318i,-0.1815+0.1318i,0.6351-0.4615i;
            0.1392-0.1759i,0.4873-0.6156i,-0.4873+0.6156i,-0.1392+0.1759i;
            0,0,0,0;
            0,0,0,0];
         
PAR.CB(:,:,4)=[  0,0,0,0;
            0,0,0,0;
            0.7851,-0.2243,0.2243,-0.7851;
            -0.0055-0.2242i,-0.0193-0.7848i,0.0193+0.7848i,0.0055+0.2242i];
 
PAR.CB(:,:,5)=[  -0.0055-0.2242i,-0.0193-0.7848i,0.0193+0.7848i,0.0055+0.2242i;
            0,0,0,0;
            0,0,0,0;
             -0.6351+0.4615i,0.1815-0.1318i,-0.1815+0.1318i,0.6351-0.4615i;];
  
PAR.CB(:,:,6)=[  0,0,0,0;
            0.7851,-0.2243,0.2243,-0.7851;
            0.1392-0.1759i,0.4873-0.6156i,-0.4873+0.6156i,-0.1392+0.1759i;
            0,0,0,0]; 
%}
  
%MC1=[0.3640 1.0500 -0.3640 -1.0500]; %2
%MC2=[-0.6068 0.6299 0.6068 -0.6299];
%MC1=[-0.9484 -0.0222 0.9484 0.0222];%2 较差
%MC2=[-0.0201 1.0486 0.0201 -1.0486];
MC1=[-1.3646 1.3646 -0.2206 0.2206];   %3
MC2=[0.3013 -0.3013 -0.9992 0.9992];
%MC1=[0.7878 -0.7370 0.7370 -0.7878];  %3 bad
%MC2=[1.2765 0.4548 -0.4548 -1.2765];
theta1=1.2;theta2=2.2;theta3=1.1;theta4=2.2;
a1=MC1;a2=MC2*(cos(theta4)+1i*sin(theta4));
b1=MC1;b2=MC2*(cos(theta4)+1i*sin(theta4));
c1=MC1*(cos(theta1)+1i*sin(theta1));c2=MC2*(cos(theta2)+1i*sin(theta2));
d1=MC1*(cos(theta3)+1i*sin(theta3));d2=MC2;
e1=MC2*(cos(theta2)+1i*sin(theta2));e2=MC1*(cos(theta3)+1i*sin(theta3));
f1=MC1*(cos(theta1)+1i*sin(theta1));f2=MC2;

PAR.CB(:,:,1)=[  0,0,0,0;
            a1;
            0,0,0,0;
            a2];
 
PAR.CB(:,:,2)=[  b1;
            0,0,0,0;
            b2;
            0,0,0,0];
         
PAR.CB(:,:,3)=[  c1;
            c2;
            0,0,0,0;
            0,0,0,0];
         
PAR.CB(:,:,4)=[  0,0,0,0;
            0,0,0,0;
            d1;
           d2];
 
PAR.CB(:,:,5)=[ e1;
            0,0,0,0;
            0,0,0,0;
             e2];
  
PAR.CB(:,:,6)=[  0,0,0,0;
            f1;
            f2;
            0,0,0,0]; 



data_source = ceil(rand(PAR.VN,PAR.Data_length)*4);    %6X4096 每列为个用户发送码字 共4096次


PRE_o = scmaEncode(data_source, PAR.CB); %4X4096 每列为6个用户发送码字在4个资源上的叠加，每个数字表示使用第几个码字


err_sym_sum(1, length(PAR.EbNo)) = 0;
SER(1, length(PAR.EbNo)) = 0;
for n = 1:20%length(PAR.EbNo)

    SNR=PAR.EbNo(n);
    Pnoise=0.667/4/(2^(SNR/10));
    %Pnoise=1/4/(2^(SNR/10));

    PRE=channel(PRE_o,Pnoise);
    total_error=0;
    for data=1:PAR.Data_length
         [err_sum,decode] = Mpa(PAR.CB, PRE(:,data), data_source(:,data),PAR,Pnoise);
         total_error = total_error+err_sum;
    end
    err_sym_sum(1, n) = total_error;
    fprintf('EbNo_dB = %d , err_sum = %d \n',PAR.EbNo(n),total_error);
   SER(1,n)=err_sym_sum(1, n)/PAR.VN/PAR.Data_length;
end

