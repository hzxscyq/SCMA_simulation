function [err_sum,decode] = Mpa(CB, PRE, data_source,PAR,Pnoise,wid)
iterNum=5;
J=size(CB,3); %J���û�
K=size(CB,1); %K����Դ�ڵ�
M=size(CB,2);   %M������
discPrec=0.1; %�������
noiseWid=5; %��������
lvalid=[2,3,5;
        1,3,6;
        2,4,6;
        1,4,5];%ÿ����Դ�ڵ����ӵĲ�ڵ�
rvalid=[2,4;
        1,3;
        1,2;
        3,4;            
        1,4;
        2,3];
decode=zeros(size(data_source));
decision=ones(J,iterNum);
V=ones(J,PAR.d_v,M)/M;
y=zeros(K,1);
discNum=pow2(nextpow2((2*(PAR.d_f-1)*wid+noiseWid)/discPrec+1));  %
noiseTmp=-noiseWid:discPrec:noiseWid;
discNoise=bsxfun(@plus,noiseTmp.^2,(noiseTmp.^2)');%generate x^2+y^2 for 4*4 square in x,y coordinate
discNoise=exp(-discNoise/(2*Pnoise))/(pi*2*Pnoise); %NXN��������

for iter=1:iterNum
    U=zeros(K,PAR.d_f,M);
    for k=1:K
        y(k)=PRE(k);
        for idf=1:PAR.d_f
            iUser=lvalid(k,idf); %��ǰ����Ĳ�ڵ�iUser
            layer=lvalid(k,[1:idf-1 idf+1:PAR.d_f]);
            probDisc=ones(discNum); %��ǰ��ڵ��
            for c=1:PAR.d_f-1
                jj=layer(c); %��iUser�������ڵ�
                temp=zeros(discNum);%��ڵ�j����ɢ�����ܶ�
                for m=1:M %��ڵ�j��M�����ָ���
                    realDisc=round((real(CB(k,m,jj))+wid)/discPrec)+1;
                    imagDisc=round((imag(CB(k,m,jj))+wid)/discPrec)+1;
                    temp(realDisc,imagDisc)=temp(realDisc,imagDisc)+V(jj,rvalid(jj,:)==k,m);
                end
                probDisc=probDisc.*fft2(temp); %df-1�������ܶȾ��
            end
            probDisc=ifft2(probDisc.*fft2(discNoise,discNum,discNum)); %��������
            for m=1:M  %����m�����ֵĸ��� ����y_k-x_kjm
                realDisc=round((real(y(k)-CB(k,m,iUser))+wid*(PAR.d_f-1)+noiseWid)/discPrec)+1;
                imagDisc=round((imag(y(k)-CB(k,m,iUser))+wid*(PAR.d_f-1)+noiseWid)/discPrec)+1;
                U(k,idf,m)=probDisc(realDisc,imagDisc);
            end
        end
    end
    
    %����V
    V=ones(J,PAR.d_v,M);
    for j=1:J  %������ڵ�
        for v=1:PAR.d_v  %������Դ�ڵ�
            resource=rvalid(j,[1:v-1 v+1:PAR.d_v]);
            for m=1:M
                for c=1:PAR.d_v-1
                    V(j,v,m)=V(j,v,m)*U(resource(c),lvalid(resource(c),:)==j,m);
                end
            end
            V(j,v,:)=V(j,v,:)/sum(V(j,v,:));
        end
    end
        %�����о�
    result=ones(J,M);
    for j=1:J
        for m=1:M
            for v=1:PAR.d_v
                result(j,m)=result(j,m)*U(rvalid(j,v),lvalid(rvalid(j,v),:)==j,m);
            end
        end
        [~,ind]=max(result(j,:));
        decision(j,iter)=ind;
    end
    
end
decode=decision(:,iterNum);
err=decode~=data_source;
err_sum=sum(sum(err));
end