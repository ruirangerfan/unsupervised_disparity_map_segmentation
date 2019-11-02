%%%%%%%%%%%%---TITS source code---%%%%%%%%%%%%
% This code was created by Dr. Rui Ranger Fan 
% email: rui.fan@ieee.org
% web:   ruirangerfan.com
%% Road Damage Detection Based on Unsupervised Disparity Map Segmentation

clc;
clear all; 
close all; 

load('data.mat')

%% disparity transformation
tic;
[vmax,umax]=size(disp);
v_map=[1:vmax]'*ones(1,umax);
u_map= ones(vmax,1)*[1:umax];

v=reshape(v_map,[vmax*umax,1]);
u=reshape(u_map,[vmax*umax,1]);
d=reshape(disp,[vmax*umax,1]);
    
u(d==0)=[];
v(d==0)=[];
d(d==0)=[];
    
Su=sum(u);
Sv=sum(v);
Sd=sum(d);
Su2=sum(u.^2);
Sv2=sum(v.^2);
Sdu=sum(u.*d);
Sdv=sum(v.*d);
Suv=sum(u.*v);
n=length(u);
 
beta0=(Sd^2*(Sv2+Su2)-2*Sd*(Sv*Sdv+Su*Sdu)+n*(Sdv^2+Sdu^2))/2;
beta1=(Sd^2*(Sv2-Su2)+2*Sd*(Su*Sdu-Sv*Sdv)+n*(Sdv^2-Sdu^2))/2;
beta2=-Sd^2*Suv+Sd*(Sv*Sdu+Su*Sdv)-n*Sdv*Sdu;
gamma0=(n*Sv2+n*Su2-Sv^2-Su^2)/2;
gamma1=(n*Sv2-n*Su2-Sv^2+Su^2)/2;
gamma2=Sv*Su-n*Suv;
clear Su Sv Sd Su2 Sv2 Sdu Sdv Suv Suv vmax umax

A=(beta1*gamma0-beta0*gamma1);
B=(beta0*gamma2-beta2*gamma0);
C=(beta1*gamma2-beta2*gamma1);
    
delta=A^2+B^2-C^2;
tmp1=(A+sqrt(delta))/(B-C);
tmp2=(A-sqrt(delta))/(B-C);
theta1=atan(tmp1);
theta2=atan(tmp2);
clear A B C beta0 beta1 beta2 gamma0 gamma1 gamma2 tmp1 tmp2 delta

t1=v*cos(theta1)-u*sin(theta1);
t2=v*cos(theta2)-u*sin(theta2);

T1=[ones(n,1),t1];
T2=[ones(n,1),t2];
    
f1=d'*T1*inv(T1'*T1)*T1'*d;
f2=d'*T2*inv(T2'*T2)*T2'*d;

if f1<f2
    theta=theta2;
else
    theta=theta1;
end
clear t1 t2 T1 T2 f1 f2 theta1 theta2

t=v*cos(theta)-u*sin(theta);
T=[ones(n,1),t];
a=inv(T'*T)*T'*d;

t_map=v_map*cos(theta)-u_map*sin(theta);
disp1=disp-(a(1)+a(2)*t_map)+30;
disp1(disp==0)=0;
clear t_map v_map u_map theta a u v d t n T 
% disp1 is the transformed disparity map
toc;

figure;
ax = subplot(2,1,1); imshow(disp, [],'Colormap',jet(4096));
title('original disparity map');
ax = subplot(2,1,2); imshow(disp1,[],'Colormap',jet(4096));
title('transformed disparity map');