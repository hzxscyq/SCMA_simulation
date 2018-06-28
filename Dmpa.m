function [err_sum,decode] = Mpa(CB, PRE, data_source,PAR,Pnoise,wid)
iterNum=5;
J=size(CB,3); %J个用户
K=size(CB,1); %K个资源节点
M=size(CB,2);   %M个码字
discPrec=0.1; %采样间隔
noiseWid=5; %噪声区间
lvalid=[2,3,5;
        1,3,6;
        2,4,6;
        1,4,5];%每个资源节点连接的层节点
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
discNoise=exp(-discNoise/(2*Pnoise))/(pi*2*Pnoise); %NXN噪声矩阵

for iter=1:iterNum
    U=zeros(K,PAR.d_f,M);
    for k=1:K
        y(k)=PRE(k);
        for idf=1:PAR.d_f
            iUser=lvalid(k,idf); %当前计算的层节点iUser
            layer=lvalid(k,[1:idf-1 idf+1:PAR.d_f]);
            probDisc=ones(discNum); %当前层节点的
            for c=1:PAR.d_f-1
                jj=layer(c); %除iUser外其余层节点
                temp=zeros(discNum);%层节点j的离散概率密度
                for m=1:M %层节点j的M个码字概率
                    realDisc=round((real(CB(k,m,jj))+wid)/discPrec)+1;
                    imagDisc=round((imag(CB(k,m,jj))+wid)/discPrec)+1;
                    temp(realDisc,imagDisc)=temp(realDisc,imagDisc)+V(jj,rvalid(jj,:)==k,m);
                end
                probDisc=probDisc.*fft2(temp); %df-1个概率密度卷积
            end
            probDisc=ifft2(probDisc.*fft2(discNoise,discNum,discNum)); %乘上噪声
            for m=1:M  %计算m个码字的概率 带入y_k-x_kjm
                realDisc=round((real(y(k)-CB(k,m,iUser))+wid*(PAR.d_f-1)+noiseWid)/discPrec)+1;
                imagDisc=round((imag(y(k)-CB(k,m,iUser))+wid*(PAR.d_f-1)+noiseWid)/discPrec)+1;
                U(k,idf,m)=probDisc(realDisc,imagDisc);
            end
        end
    end
    
    %计算V
    V=ones(J,PAR.d_v,M);
    for j=1:J  %遍历层节点
        for v=1:PAR.d_v  %遍历资源节点
            resource=rvalid(j,[1:v-1 v+1:PAR.d_v]);
            for m=1:M
                for c=1:PAR.d_v-1
                    V(j,v,m)=V(j,v,m)*U(resource(c),lvalid(resource(c),:)==j,m);
                end
            end
            V(j,v,:)=V(j,v,:)/sum(V(j,v,:));
        end
    end
        %进行判决
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