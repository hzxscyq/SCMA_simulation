%% SCMA simple Transceiverchain.m
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
PAR.Data_length=4096;    % Number of total data
PAR.N_iter= 5;  % Number of iterations in decoding
%PAR.EbNo = 10;
PAR.EbNo = [1:0.1:20];
PAR.CB=zeros(PAR.M,PAR.FN,PAR.VN);      %Codebooks
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
        
%% AWGN Channel
hChan = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (Eb/No)',...
  'SignalPower', 1, 'SamplesPerSymbol', 1,'BitsPerSymbol',2);

%% Initialization
data_source = ceil(rand(PAR.VN,PAR.Data_length)*4);    %6X4096 每列为个用户发送码字 共4096次
data_source_2 = four2two(data_source);


PRE_o = scmaEncode(data_source, PAR.CB); %4X4096 每列为6个用户发送码字在4个资源上的叠加，每个数字表示使用第几个码字

%% addition AWGN and deconding
%PRE = step(hChan, PRE_o);
err_sym_sum(1, length(PAR.EbNo)) = 0;
SER(1, length(PAR.EbNo)) = 0;
for n = 1:length(PAR.EbNo)
    reset(hChan)    
        hChan.EbNo = PAR.EbNo(n);
    PRE = step(hChan, PRE_o);   %添加噪声后的接收信号
    
    % MPA detect
    [err_sum,demodSignal] = scmaDeML(PAR.CB, PRE, data_source);

    
    err_sym_sum(1, n) = err_sum;
    fprintf('EbNo_dB = %d , err_sum = %d \n',PAR.EbNo(n),err_sum);
    SER(1,n)=err_sym_sum(1, n)/PAR.VN/PAR.Data_length;
end

semilogy(PAR.EbNo, SER,'b+-');grid;xlabel('Eb/No (dB)');ylabel('BER');
