clc
clear all
close all

% Simulations settings
dist=0;         % 1 with payload / disturbances; 0 without

T=5;

n=4;

if n==3     % support file
[t,y] = ode23(@(t,y)Topics_in_control_ode_n3(t,y,dist),[0 T],[0 0 0 0 0 0]);
end

if n==4     % coursework assignment
[t,y] = ode23(@(t,y)Topics_in_control_ode_n4_solutions(t,y,dist),[0 T],[0 0 0 0 0 0 0 0]);   
end

xd=pi/6;
ref=repmat(xd,size(t));  


%% Figures

% Position plot
if n==3
    
disp('total error [rad]');
disp(y(end,1)+y(end,3)+y(end,5)-3*ref(end)); 
    
h1=figure(1);
plot(t,y(:,1),'b','LineWidth',2);hold on;plot(t,y(:,3),'--r','LineWidth',2);plot(t,y(:,5),':g','LineWidth',2);
plot(t,ref,':k','LineWidth',2);
hold off;
xlabel('time [s]');
ylabel('Joint angles [rad]');    
legend('q1','q2','q3');
set(h1, 'Position', [100, 100, 350, 300]);
axis([0 T 0 0.6]);
end

if n==4
    
disp('total error [rad]');
disp(y(end,1)+y(end,3)+y(end,5)+y(end,7)-4*ref(end)); 

h1=figure(1);
plot(t,y(:,1),'b','LineWidth',2);hold on;plot(t,y(:,3),'--r','LineWidth',2);plot(t,y(:,5),':g','LineWidth',2);plot(t,y(:,7),'-.m','LineWidth',2);
plot(t,ref,':k','LineWidth',2);
hold off;
xlabel('time [s]');
ylabel('Joint angles [rad]');    
legend('q1','q2','q3','q4');
set(h1, 'Position', [100, 100, 350, 300]);
axis([0 T 0 0.6]);
end

