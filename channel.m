function [ PRE] = channel( PRE_o,nVar )
%c[J][K] is the K-length codeword from the J-th layer
%nVar is the noise variance in every dimensions
%y[K] is the K-length received codeword
%h[J][K] is the channel vectors for the J layers 

K=size(PRE_o,1);
length = size(PRE_o,2);
PRE=zeros(K,length);


for len=1:length
    for i=1:K
        noise=(randn+1i*randn)*sqrt(nVar);
        PRE(i,len)=PRE_o(i,len)+noise;
    end
end




end

