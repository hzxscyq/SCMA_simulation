clc;clear;
%MC1=[-0.8621 0.8621 0.1429 -0.1429];
%MC2=[0.0825 -0.0825 1.4932 -1.4932];
MC1=[-1.3646 1.3646 -0.2206 0.2206];
MC2=[0.3013 -0.3013 -0.9992 0.9992];%31
%MC1=[-0.9484 -0.0222 0.9484 0.0222];
%MC2=[-0.0201 1.0486 0.0201 -1.0486]; %2 bad
%MC1=[0.3640 1.0500 -0.3640 -1.0500]; %2
%MC2=[-0.6068 0.6299 0.6068 -0.6299];
%MC1=[0.7878 -0.7370 0.7370 -0.7878];  %3 bad
%MC2=[1.2765 0.4548 -0.4548 -1.2765];
C1=MC2;
distance=0;
theta_temp=0;
theta2_temp=0;
for theta=0:0.1:pi
    for theta2=0:0.1:pi
        C2=MC1*(cos(theta)+1i*sin(theta));
        C3=MC2*(cos(theta2)+1i*sin(theta2));
        dis=distance_cal(C1,C2,C3);
        if dis>distance
            distance=dis;
            theta_temp=theta;
            theta2_temp=theta2;
        end
    end
end
