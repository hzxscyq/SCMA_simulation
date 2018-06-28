clc;clear;
N=2;M=4;
IN=eye(N);
IM=eye(M);
Ei=cell(4,1);
E=cell(6,1);
e=[1 0 0 0 0 0 0 0;
   0 1 0 0 0 0 0 0;
   0 0 -1 0 0 0 0 0;
   0 0 0 -1 0 0 0 0;
   0 0 0 0 0 0 0 0;
   0 0 0 0 0 0 0 0;
   0 0 0 0 0 0 0 0;
   0 0 0 0 0 0 0 0];
m=1;
D=1.3;
for i=1:M
    Ei{i}=kron(IM(:,i)',IN);
end

for i=1:M
    for j=i+1:M
        E{m}=Ei{i}'*Ei{i}-Ei{i}'*Ei{j}-Ei{j}'*Ei{i}+Ei{j}'*Ei{j};
        m=m+1;
    end
end
x=randn(M*N,1);

flag=test(x,E,D);
while flag==1
    x=randn(M*N,1);
    flag=test(x,E,D);
end
%获得初始数据点x0
count=1;
cvx_begin quiet
    variable x1(M*N)
    minimize(norm(x1,2))
    subject to
        2*x'*E{1}*x1-x'*E{1}*x-D>=0
        2*x'*E{2}*x1-x'*E{2}*x-D>=0
        2*x'*E{3}*x1-x'*E{3}*x-D>=0
        2*x'*E{4}*x1-x'*E{4}*x-D>=0
        2*x'*E{5}*x1-x'*E{5}*x-D>=0
        2*x'*E{6}*x1-x'*E{6}*x-D>=0 
       % x1(1)^2+x1(2)^2-x1(3)^2-x1(4)^2 == 0
        
        %{
        2*x'*E{7}*x1-x'*E{7}*x-D>=0
        2*x'*E{8}*x1-x'*E{8}*x-D>=0
        2*x'*E{9}*x1-x'*E{9}*x-D>=0
        2*x'*E{10}*x1-x'*E{10}*x-D>=0
        2*x'*E{11}*x1-x'*E{11}*x-D>=0
        2*x'*E{12}*x1-x'*E{12}*x-D>=0
        2*x'*E{13}*x1-x'*E{13}*x-D>=0
        2*x'*E{14}*x1-x'*E{14}*x-D>=0
        2*x'*E{15}*x1-x'*E{15}*x-D>=0
        %}
cvx_end

distance(count)=norm(x1,2);
count=count+1;
while(norm(x1-x,2)>=0.01)
    x=x1;
    cvx_begin quiet
        variable x1(M*N)
        minimize(norm(x1,2))
        subject to
            2*x'*E{1}*x1-x'*E{1}*x-D>=0
            2*x'*E{2}*x1-x'*E{2}*x-D>=0
            2*x'*E{3}*x1-x'*E{3}*x-D>=0
            2*x'*E{4}*x1-x'*E{4}*x-D>=0
            2*x'*E{5}*x1-x'*E{5}*x-D>=0
            2*x'*E{6}*x1-x'*E{6}*x-D>=0
           
    cvx_end
    norm(x1-x,2)
    distance(count)=norm(x1,2);
    count=count+1;
end
plot(distance,'bs-');