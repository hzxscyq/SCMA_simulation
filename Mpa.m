function [err_sum,decode] = Mpa(CB, PRE, data_source,PAR,Pnoise)
J=size(CB,3); %J���û�
K=size(CB,1); %K����Դ�ڵ�
M=size(CB,2);   %M������
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
    
iterNum=4;%ѭ��5��
decode=zeros(size(data_source));
decision=ones(J,iterNum);
V=ones(J,PAR.d_v,M)/M;
y=zeros(K,1);

for iter=1:iterNum
    U=zeros(K,PAR.d_f,M);
    for k=1:K %������Դ�ڵ�
        y(k)=PRE(k);
        for idf=1:PAR.d_f %�����뵱ǰ��Դ�ڵ����������в�ڵ�
            iUser=lvalid(k,idf);%��ǰ�����ڵ�
            layer=lvalid(k,[1:idf-1 idf+1:PAR.d_f]); %����ǰ��ڵ�������df-1����ڵ����
            for m=1:M  %������ǰ��ڵ����������
                for n=0:M^(PAR.d_f-1)-1 %���Ľ��Ʊ�ʾ 15-33 ��000��M-1 M-1 M-1������
                    mset=dec2base(n,M,PAR.d_f-1);%ʮ������nתΪM���ƣ�����df-1������
                    tmp=y(k)-CB(k,m,iUser); %y_k-x_kjm
                    for c=1:PAR.d_f-1
                        tmp=tmp-CB(k,base2dec(mset(c),M)+1,layer(c)); %y_k-x_kjm-c_i
                    end
                    tmp=exp(-(tmp*tmp')/(2*Pnoise))/(2*pi*Pnoise);
                    for c=1:PAR.d_f-1 %����V
                        tmp=tmp*V(layer(c),rvalid(layer(c),:)==k,base2dec(mset(c),M)+1);
                    end
                    U(k,idf,m)=U(k,idf,m)+tmp;
                end
            end

        end
    end
    
    %����V
    V=ones(J,PAR.d_v,M);
    for j=1:J %�������в�ڵ�
        for v=1:PAR.d_v %����������Դ�ڵ�
            resource=rvalid(j,[1:v-1 v+1:PAR.d_v]);
            for m=1:M
                for vm=1:PAR.d_v-1
                    V(j,v,m)=V(j,v,m)*U(resource(vm),lvalid(resource(vm),:)==j,m);
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

