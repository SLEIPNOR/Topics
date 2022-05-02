function dydt = Topics_in_control_ode_n4_solutions(t,y,dist)

dydt = zeros(8,1);    
x1=y(1);
v1=y(2);
x2=y(3);
v2=y(4);
x3=y(5);
v3=y(6);
x4=y(7);
v4=y(8);
% model parameters

l=0.025;          % in [m]


m=4;        % in [g]


k=5;                % in [mNm/rad]
b=0.01;             % physical damping [Nmm*s/rad]
xd=pi/6;

% joint angles defined as relative rotations from previous link (as in potential energy)
M =  [(l^2*m*(4*cos(x2 + x3 + x4) + 12*cos(x2 + x3) + 4*cos(x3 + x4) + 20*cos(x2) + 12*cos(x3) + 4*cos(x4) + 28))/4, (l^2*m*(2*cos(x2 + x3 + x4) + 6*cos(x2 + x3) + 4*cos(x3 + x4) + 10*cos(x2) + 12*cos(x3) + 4*cos(x4) + 15))/4, (l^2*m*(cos(x2 + x3 + x4) + 3*cos(x2 + x3) + cos(x3 + x4) + 3*cos(x3) + 2*cos(x4) + 3))/2, (l^2*m*(2*cos(x2 + x3 + x4) + 2*cos(x3 + x4) + 2*cos(x4) + 1))/4;(l^2*m*(2*cos(x2 + x3 + x4) + 6*cos(x2 + x3) + 4*cos(x3 + x4) + 10*cos(x2) + 12*cos(x3) + 4*cos(x4) + 15))/4, (l^2*m*(4*cos(x3 + x4) + 12*cos(x3) + 4*cos(x4) + 15))/4, (l^2*m*(cos(x3 + x4) + 3*cos(x3) + 2*cos(x4) + 3))/2, (l^2*m*(2*cos(x3 + x4) + 2*cos(x4) + 1))/4;(l^2*m*(cos(x2 + x3 + x4) + 3*cos(x2 + x3) + cos(x3 + x4) + 3*cos(x3) + 2*cos(x4) + 3))/2, (l^2*m*(cos(x3 + x4) + 3*cos(x3) + 2*cos(x4) + 3))/2, (l^2*m*(4*cos(x4) + 6))/4, (l^2*m*(2*cos(x4) + 1))/4;(l^2*m*(2*cos(x2 + x3 + x4) + 2*cos(x3 + x4) + 2*cos(x4) + 1))/4, (l^2*m*(2*cos(x3 + x4) + 2*cos(x4) + 1))/4, (l^2*m*(2*cos(x4) + 1))/4, (l^2*m)/4];
dMt =  [- (v4*l^2*m*(4*sin(x2 + x3 + x4) + 4*sin(x3 + x4) + 4*sin(x4)))/4 - (v2*l^2*m*(4*sin(x2 + x3 + x4) + 12*sin(x2 + x3) + 20*sin(x2)))/4 - (v3*l^2*m*(4*sin(x2 + x3 + x4) + 12*sin(x2 + x3) + 4*sin(x3 + x4) + 12*sin(x3)))/4, - (v4*l^2*m*(2*sin(x2 + x3 + x4) + 4*sin(x3 + x4) + 4*sin(x4)))/4 - (v2*l^2*m*(2*sin(x2 + x3 + x4) + 6*sin(x2 + x3) + 10*sin(x2)))/4 - (v3*l^2*m*(2*sin(x2 + x3 + x4) + 6*sin(x2 + x3) + 4*sin(x3 + x4) + 12*sin(x3)))/4, - (v4*l^2*m*(sin(x2 + x3 + x4) + sin(x3 + x4) + 2*sin(x4)))/2 - (v3*l^2*m*(sin(x2 + x3 + x4) + 3*sin(x2 + x3) + sin(x3 + x4) + 3*sin(x3)))/2 - (v2*l^2*m*(sin(x2 + x3 + x4) + 3*sin(x2 + x3)))/2, - (v2*l^2*m*sin(x2 + x3 + x4))/2 - (v4*l^2*m*(2*sin(x2 + x3 + x4) + 2*sin(x3 + x4) + 2*sin(x4)))/4 - (v3*l^2*m*(2*sin(x2 + x3 + x4) + 2*sin(x3 + x4)))/4;- (v4*l^2*m*(2*sin(x2 + x3 + x4) + 4*sin(x3 + x4) + 4*sin(x4)))/4 - (v2*l^2*m*(2*sin(x2 + x3 + x4) + 6*sin(x2 + x3) + 10*sin(x2)))/4 - (v3*l^2*m*(2*sin(x2 + x3 + x4) + 6*sin(x2 + x3) + 4*sin(x3 + x4) + 12*sin(x3)))/4, - (v4*l^2*m*(4*sin(x3 + x4) + 4*sin(x4)))/4 - (v3*l^2*m*(4*sin(x3 + x4) + 12*sin(x3)))/4, - (v3*l^2*m*(sin(x3 + x4) + 3*sin(x3)))/2 - (v4*l^2*m*(sin(x3 + x4) + 2*sin(x4)))/2, - (v3*l^2*m*sin(x3 + x4))/2 - (v4*l^2*m*(2*sin(x3 + x4) + 2*sin(x4)))/4;- (v4*l^2*m*(sin(x2 + x3 + x4) + sin(x3 + x4) + 2*sin(x4)))/2 - (v3*l^2*m*(sin(x2 + x3 + x4) + 3*sin(x2 + x3) + sin(x3 + x4) + 3*sin(x3)))/2 - (v2*l^2*m*(sin(x2 + x3 + x4) + 3*sin(x2 + x3)))/2, - (v3*l^2*m*(sin(x3 + x4) + 3*sin(x3)))/2 - (v4*l^2*m*(sin(x3 + x4) + 2*sin(x4)))/2, -v4*l^2*m*sin(x4), -(v4*l^2*m*sin(x4))/2;- (v2*l^2*m*sin(x2 + x3 + x4))/2 - (v4*l^2*m*(2*sin(x2 + x3 + x4) + 2*sin(x3 + x4) + 2*sin(x4)))/4 - (v3*l^2*m*(2*sin(x2 + x3 + x4) + 2*sin(x3 + x4)))/4, - (v3*l^2*m*sin(x3 + x4))/2 - (v4*l^2*m*(2*sin(x3 + x4) + 2*sin(x4)))/4, -(v4*l^2*m*sin(x4))/2, 0];
GrV=[k*x1;k*x2;k*x3;k*x4];
G=[1;1;1;1];


% payload / disturbances
if dist==1
    f_dist=[1;1;1;1];  
    d0=0.1;       % tip moment
else
    f_dist=[0;0;0;0];  
    d0=0;
end

% control law

Kp=0.1;
Kv=2;
Km=5;

p=M*[v1;v2;v3;v4];

alpha_d=4*xd;

N=4;
alpha=x1+x2+x3+x4;
alpha_f=v1+v2+v3+v4;
U1=k*alpha/N-Kp*Km*(alpha-alpha_d)-Kv/Km*alpha_f;

if dist==1
  U1=U1+d0;
end
% if dist==1  %% add compensator after 3s
%     if t>3 
%         U1=U1+d0;  
%     end
% end

% Equations of motion
dqH = [k*x1;k*x2 - (l^2*m*(20*v1^2*sin(x2) + 4*v1^2*sin(x2 + x3 + x4) + 12*v1^2*sin(x2 + x3) + 20*v1*v2*sin(x2) + 4*v1*v2*sin(x2 + x3 + x4) + 4*v1*v3*sin(x2 + x3 + x4) + 4*v1*v4*sin(x2 + x3 + x4) + 12*v1*v2*sin(x2 + x3) + 12*v1*v3*sin(x2 + x3)))/8;k*x3 - (l^2*m*(12*v1^2*sin(x3) + 12*v2^2*sin(x3) + 4*v1^2*sin(x2 + x3 + x4) + 12*v1^2*sin(x2 + x3) + 4*v1^2*sin(x3 + x4) + 4*v2^2*sin(x3 + x4) + 24*v1*v2*sin(x3) + 12*v1*v3*sin(x3) + 12*v2*v3*sin(x3) + 4*v1*v2*sin(x2 + x3 + x4) + 4*v1*v3*sin(x2 + x3 + x4) + 4*v1*v4*sin(x2 + x3 + x4) + 12*v1*v2*sin(x2 + x3) + 12*v1*v3*sin(x2 + x3) + 8*v1*v2*sin(x3 + x4) + 4*v1*v3*sin(x3 + x4) + 4*v1*v4*sin(x3 + x4) + 4*v2*v3*sin(x3 + x4) + 4*v2*v4*sin(x3 + x4)))/8;k*x4 - (l^2*m*(4*v1^2*sin(x4) + 4*v2^2*sin(x4) + 4*v3^2*sin(x4) + 4*v1^2*sin(x2 + x3 + x4) + 4*v1^2*sin(x3 + x4) + 4*v2^2*sin(x3 + x4) + 8*v1*v2*sin(x4) + 8*v1*v3*sin(x4) + 4*v1*v4*sin(x4) + 8*v2*v3*sin(x4) + 4*v2*v4*sin(x4) + 4*v3*v4*sin(x4) + 4*v1*v2*sin(x2 + x3 + x4) + 4*v1*v3*sin(x2 + x3 + x4) + 4*v1*v4*sin(x2 + x3 + x4) + 8*v1*v2*sin(x3 + x4) + 4*v1*v3*sin(x3 + x4) + 4*v1*v4*sin(x3 + x4) + 4*v2*v3*sin(x3 + x4) + 4*v2*v4*sin(x3 + x4)))/8];
dpH = (M^-1)*p;
dp = -dqH+G*U1-d0*f_dist-b*dpH;
dv=(M^-1)*(dp-dMt*[v1;v2;v3;v4]);

dydt(1)=y(2);
dydt(2)=dv(1);
dydt(3)=y(4);
dydt(4)=dv(2);
dydt(5)=y(6);
dydt(6)=dv(3);
dydt(7)=y(8);
dydt(8)=dv(4);
