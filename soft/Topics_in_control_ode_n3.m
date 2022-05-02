function dydt = Topics_in_control_ode_n3(t,y,dist)

dydt = zeros(6,1);    
x1=y(1);
v1=y(2);
x2=y(3);
v2=y(4);
x3=y(5);
v3=y(6);

% model parameters

l0=0.1;             % in [m]
l1=l0*1/100;
l2=l0*33/100;
l3=l0*33/100;
l4=l0*33/100;

m0=15;              % in [g]
m1=m0*1/100;
m2=m0*33/100;
m3=m0*33/100;
m4=m0*33/100;

k=5;                % in [mNm/rad]
b=0.01;             % physical damping [Nmm*s/rad]
xd=pi/6;

% joint angles defined as relative rotations from previous link (as in potential energy)
M=[(l2^2*m2)/4 + l2^2*m3 + l2^2*m4 + (l3^2*m3)/4 + l3^2*m4 + (l4^2*m4)/4 + l2*l4*m4*cos(x2 + x3) + l2*l3*m3*cos(x2) + 2*l2*l3*m4*cos(x2) + l3*l4*m4*cos(x3), (l3^2*m3)/4 + l3^2*m4 + (l4^2*m4)/4 + (l2*l4*m4*cos(x2 + x3))/2 + l3*l4*m4*cos(x3) + l2*l3*cos(x2)*(m3/2 + m4), (l4^2*m4)/4 + (l2*l4*m4*cos(x2 + x3))/2 + (l3*l4*m4*cos(x3))/2;(l3^2*m3)/4 + l3^2*m4 + (l4^2*m4)/4 + (l2*l4*m4*cos(x2 + x3))/2 + l3*l4*m4*cos(x3) + l2*l3*cos(x2)*(m3/2 + m4), (l3^2*m3)/4 + l3^2*m4 + (l4^2*m4)/4 + l3*l4*m4*cos(x3), (m4*l4^2)/4 + (l3*m4*cos(x3)*l4)/2;(l4^2*m4)/4 + (l2*l4*m4*cos(x2 + x3))/2 + (l3*l4*m4*cos(x3))/2, (m4*l4^2)/4 + (l3*m4*cos(x3)*l4)/2, (l4^2*m4)/4];
m11=-l2*l4*m4*sin(x2 + x3)*(v2+v3) - l2*l3*m3*sin(x2)*v2 -2*l2*l3*m4*sin(x2)*v2 - l3*l4*m4*sin(x3)*v3;
dMt=[m11, -l2*l4*m4*sin(x2 + x3)*(v2+v3)/2 - l3*l4*m4*sin(x3)*v3 - l2*l3*sin(x2)*(m3/2 + m4)*v2, -l2*l4*m4*sin(x2 + x3)*(v2+v3)/2 - l3*l4*m4*sin(x3)*v3/2; -l2*l4*m4*sin(x2 + x3)*(v2+v3)/2 - l3*l4*m4*sin(x3)*v3 - l2*l3*sin(x2)*(m3/2 + m4)*v2, - l3*l4*m4*sin(x3)*v3, -l3*m4*sin(x3)*l4*v3/2; -l2*l4*m4*sin(x2 + x3)*(v2+v3)/2 - l3*l4*m4*sin(x3)*v3/2, -l3*m4*sin(x3)*l4*v3/2, 0];

GrV=[k*x1;k*x2;k*x3];
G=[1;1;1];


% payload / disturbances
if dist==1
    f_dist=[1;1;1];  
    d0=0.1;       % tip moment
else
    f_dist=[0;0;0];  
    d0=0;
end

% control law

Kp=0.1;
Kv=2;
Km=5;

p=M*[v1;v2;v3];

alpha_d=3*xd;

N=3;
alpha=x1+x2+x3;
alpha_f=v1+v2+v3;
U1=k*alpha/N-Kp*Km*(alpha-alpha_d)-Kv/Km*alpha_f;

if dist==1
  U1=U1+d0;  
end

% Equations of motion
dp=-[0;- v1*(v2*((l2*l4*m4*sin(x2 + x3))/2 + l2*l3*sin(x2)*(m3/2 + m4)) + v1*(l2*l4*m4*sin(x2 + x3) + l2*l3*m3*sin(x2) + 2*l2*l3*m4*sin(x2)) + (l2*l4*m4*v3*sin(x2 + x3))/2) - v1*v2*((l2*l4*m4*sin(x2 + x3))/2 + l2*l3*sin(x2)*(m3/2 + m4)) - (l2*l4*m4*v1*v3*sin(x2 + x3))/2;- v2*(v1*((l2*l4*m4*sin(x2 + x3))/2 + l3*l4*m4*sin(x3)) + l3*l4*m4*v2*sin(x3) + (l3*l4*m4*v3*sin(x3))/2) - v3*(v1*((l2*l4*m4*sin(x2 + x3))/2 + (l3*l4*m4*sin(x3))/2) + (l3*l4*m4*v2*sin(x3))/2) - v1*(v1*(l2*l4*m4*sin(x2 + x3) + l3*l4*m4*sin(x3)) + v2*((l2*l4*m4*sin(x2 + x3))/2 + l3*l4*m4*sin(x3)) + v3*((l2*l4*m4*sin(x2 + x3))/2 + (l3*l4*m4*sin(x3))/2))]-GrV+G*U1-d0*f_dist-b*[v1;v2;v3];
iM=[(4*m3 + 8*m4 - 8*m4*cos(2*x3))/(l2^2*(m2*m3 + 2*m2*m4 + 6*m3*m4 - 2*m3^2*cos(2*x2) + 2*m3^2 - 2*m2*m4*cos(2*x3) - 4*m3*m4*cos(2*x2) - 4*m3*m4*cos(2*x3) + 2*m3*m4*cos(2*x2 + 2*x3))), -(4*l3*m3 + 8*l3*m4 + 8*l2*m3*cos(x2) + 8*l2*m4*cos(x2) - 8*l2*m4*cos(x2 + 2*x3) - 8*l3*m4*cos(2*x3))/(l2^2*l3*(m2*m3 + 2*m2*m4 + 10*m3*m4 - 4*m3^2*cos(x2)^2 + 4*m3^2 - 2*m2*m4*cos(2*x3) - 4*m3*m4*cos(2*x3) - 8*m3*m4*cos(x2)^2 + 2*m3*m4*cos(2*x2 + 2*x3))), (8*l4*m3*cos(x2) + 8*l4*m4*cos(x2) + 8*l3*m3*cos(x2 - x3) + 16*l3*m4*cos(x2 - x3) - 8*l4*m4*cos(x2 + 2*x3) - 16*l3*m4*cos(x2 + x3))/(l2*l3*l4*(m2*m3 + 2*m2*m4 + 10*m3*m4 - 4*m3^2*cos(x2)^2 + 4*m3^2 - 2*m2*m4*cos(2*x3) - 4*m3*m4*cos(2*x3) - 8*m3*m4*cos(x2)^2 + 2*m3*m4*cos(2*x2 + 2*x3)));-(4*l3*m3 + 8*l3*m4 + 8*l2*m3*cos(x2) + 8*l2*m4*cos(x2) - 8*l2*m4*cos(x2 + 2*x3) - 8*l3*m4*cos(2*x3))/(l2^2*l3*(m2*m3 + 2*m2*m4 + 10*m3*m4 - 4*m3^2*cos(x2)^2 + 4*m3^2 - 2*m2*m4*cos(2*x3) - 4*m3*m4*cos(2*x3) - 8*m3*m4*cos(x2)^2 + 2*m3*m4*cos(2*x2 + 2*x3))), (4*l2^2*m2 + 16*l2^2*m3 + 8*l2^2*m4 + 4*l3^2*m3 + 8*l3^2*m4 - 8*l3^2*m4*cos(2*x3) - 8*l2^2*m4*cos(2*x2 + 2*x3) + 16*l2*l3*m3*cos(x2) + 16*l2*l3*m4*cos(x2) - 16*l2*l3*m4*cos(x2 + 2*x3))/(l2^2*l3^2*(m2*m3 + 2*m2*m4 + 10*m3*m4 - 4*m3^2*cos(x2)^2 + 4*m3^2 - 2*m2*m4*cos(2*x3) - 4*m3*m4*cos(2*x3) - 8*m3*m4*cos(x2)^2 + 2*m3*m4*cos(2*x2 + 2*x3))), -(8*l3^2*m3*cos(x2 - x3) - 16*l3^2*m4*cos(x2 + x3) + 16*l3^2*m4*cos(x2 - x3) + 4*l2*l4*m2 + 16*l2*l4*m3 + 8*l2*l4*m4 - 8*l2*l4*m4*cos(2*x2 + 2*x3) + 8*l2*l3*m2*cos(x3) + 24*l2*l3*m3*cos(x3) + 16*l2*l3*m4*cos(x3) + 8*l3*l4*m3*cos(x2) + 8*l3*l4*m4*cos(x2) - 8*l2*l3*m3*cos(2*x2 + x3) - 16*l2*l3*m4*cos(2*x2 + x3) - 8*l3*l4*m4*cos(x2 + 2*x3))/(l2*l3^2*l4*(m2*m3 + 4*m2*m4 + 14*m3*m4 - 4*m3^2*cos(x2)^2 + 4*m3^2 - 4*m2*m4*cos(x3)^2 - 8*m3*m4*cos(x2)^2 - 8*m3*m4*cos(x3)^2 + 2*m3*m4*cos(2*x2 + 2*x3)));(8*l4*m3*cos(x2) + 8*l4*m4*cos(x2) + 8*l3*m3*cos(x2 - x3) + 16*l3*m4*cos(x2 - x3) - 8*l4*m4*cos(x2 + 2*x3) - 16*l3*m4*cos(x2 + x3))/(l2*l3*l4*(m2*m3 + 2*m2*m4 + 10*m3*m4 - 4*m3^2*cos(x2)^2 + 4*m3^2 - 2*m2*m4*cos(2*x3) - 4*m3*m4*cos(2*x3) - 8*m3*m4*cos(x2)^2 + 2*m3*m4*cos(2*x2 + 2*x3))), -(8*l3^2*m3*cos(x2 - x3) - 16*l3^2*m4*cos(x2 + x3) + 16*l3^2*m4*cos(x2 - x3) + 4*l2*l4*m2 + 16*l2*l4*m3 + 8*l2*l4*m4 - 8*l2*l4*m4*cos(2*x2 + 2*x3) + 8*l2*l3*m2*cos(x3) + 24*l2*l3*m3*cos(x3) + 16*l2*l3*m4*cos(x3) + 8*l3*l4*m3*cos(x2) + 8*l3*l4*m4*cos(x2) - 8*l2*l3*m3*cos(2*x2 + x3) - 16*l2*l3*m4*cos(2*x2 + x3) - 8*l3*l4*m4*cos(x2 + 2*x3))/(l2*l3^2*l4*(m2*m3 + 4*m2*m4 + 14*m3*m4 - 4*m3^2*cos(x2)^2 + 4*m3^2 - 4*m2*m4*cos(x3)^2 - 8*m3*m4*cos(x2)^2 - 8*m3*m4*cos(x3)^2 + 2*m3*m4*cos(2*x2 + 2*x3))), (8*l3^2*m3^2 + 32*l3^2*m4^2 + 8*l4^2*m4^2 + 4*l3^2*m2*m3 + 16*l3^2*m2*m4 + 48*l3^2*m3*m4 + 4*l4^2*m2*m4 + 16*l4^2*m3*m4 - 8*l3^2*m3^2*cos(2*x2) - 32*l3^2*m4^2*cos(2*x2) - 8*l4^2*m4^2*cos(2*x2 + 2*x3) - 32*l3*l4*m4^2*cos(2*x2 + x3) - 32*l3^2*m3*m4*cos(2*x2) + 32*l3*l4*m4^2*cos(x3) - 16*l3*l4*m3*m4*cos(2*x2 + x3) + 16*l3*l4*m2*m4*cos(x3) + 48*l3*l4*m3*m4*cos(x3))/(l3^2*l4^2*m4*(m2*m3 + 4*m2*m4 + 10*m3*m4 - 2*m3^2*cos(2*x2) + 2*m3^2 - 4*m3*m4*cos(2*x2) - 4*m2*m4*cos(x3)^2 - 8*m3*m4*cos(x3)^2 + 2*m3*m4*cos(2*x2 + 2*x3)))];
dv=iM*(dp-dMt*[v1;v2;v3]);

dydt(1)=y(2);
dydt(2)=dv(1);
dydt(3)=y(4);
dydt(4)=dv(2);
dydt(5)=y(6);
dydt(6)=dv(3);
