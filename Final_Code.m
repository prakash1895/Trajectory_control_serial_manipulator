function []= Throw_arm_2()
clc
clear all;
close all;
set(groot,'DefaultFigureRenderer','painters');
format long ;

% Manipulator Parameters
I1=0.125;  I2 =0.25 ; m1=3; m2=3; h=0.5; l=0.5; g=9.81; m=0;

t0 = 0;
tf = 5;

x_des = 1;
y_des = 1;

d = sqrt((x_des^2 + y_des^2));
vf = (d)/sqrt((h+l)/4.9); 

wf = vf/l; 

q1_des = atan2(-y_des,x_des);

global torque omega
omega = [];
torque = [];
wf_des = [] ;
theta_des = [] ;

x0 = [0,0,0,0]; % Initial Condition - Format:[theta1,theta2,dtheta1,dtheta2]

tf = 3;

index = 1 ;
[T1,X1] = ode45(@(t,x)planarArmODE(t,x,index),[0 tf],x0);
tf = 0:0.001:3;

x1 = [X1(end,3), X1(end,4)] ; 
options=odeset('Reltol',0.001,'Stats','on','MaxStep',0.001);
[T2,X2] = ode45(@(t,x)throwArmODE(t,x),tf,x1,options);

T = [T1;3+T2];

X3 = [] ;

X3(:,2) = X2(:,1) ;
X3(:,4) = X2(:,2) ;
X3(:,1) = X1(end,1) ;
X3(:,3) = X1(end,3) ;

X = [X1;X3];
pos = forward_kinematics(T,X);

for i=1:size(T)
    wf_des = [wf_des, wf];
    theta_des = [theta_des, pi/2];
end
%% Stick Figure Animation

 n = size(T)/500;
 n = (fix(n(1)));

figure('Name','Stick Figure');
axes;
for i = 1:n:size(T)
    plot3([0 pos(i,1) pos(i,4)], [0 pos(i,2) pos(i,5)], [0 pos(i,3) pos(i,6)], '-o', 'color', [0 0.4 0.7], 'LineWidth', 3, 'MarkerSize', 5, 'MarkerFaceColor', 'r');
    title(sprintf('Planar Arm Computed Torque control \n Time: %0.2fsec \n',T(i))) ;
    xlabel('X Axis');
    ylabel('Y Axis');
    zlabel('Z Axis');
    box on;
    grid on;
    axis equal;
    axis([-1.5 1.5 -1.5 1.5 0.0 1.5]);
    drawnow();
end

%% Plot Data

figure('Name','Input Torque under Compliance Control');
plot(T, torque(1,1:size(T,1)),'r-' );
hold on
plot(T, torque(2,1:size(T,1)),'b-');
hold on
title('Control Input vs Time') ; 
xlabel('Time (sec)') ;
ylabel('Torque (Nm)') ;
legend('\tau 1','\tau 2') ;
box on;
grid on;

figure('Name','Joint velocity under Compliance Control');
plot(T, X(:,3),'r-');
hold on
plot(T, X(:,4),'b-');
hold on
plot(T, wf_des, 'g-')
title('Joint Velocity (rad/s) vs Time') ; 
xlabel('Time (sec)') ;
ylabel('Velocity (rad/s)') ;
legend('q1 dot','q2 dot', 'q_des') ;
box on;
grid on;

figure('Name','Joint angle under Compliance Control');
plot(T, X(:,1),'r-');
hold on
plot(T, X(:,2),'b-');
hold on
plot(T, theta_des, 'g-')
title('Joint Angles (rad) vs Time') ; 
xlabel('Time (sec)') ;
ylabel('Angle (rad)') ;
legend('q1','q2') ;
box on;
grid on;

%% Definging Functions
    function pos = forward_kinematics(T,X)
        [n,c] = size(T);
        for i=1:n
            pos(i,1)= 0;
            pos(i,2)= 0;
            pos(i,3)= h;
            pos(i,4)= l*cos(X(i,2))*cos(X(i,1));
            pos(i,5)= l*cos(X(i,2))*sin(X(i,1));
            pos(i,6)= h + l*sin(X(i,2));
        end
    end
     
    function dx = planarArmODE(t,x,index)
        if (index == 1)
        
            theta_d = [q1_des; 0];
            dtheta_d=[0; 0]; % Desired velocity (Derivative of theta_d)
            ddtheta_d=[0; 0];
        
        elseif (index == 2)
            
            theta_d = [q1_des; 0];
            dtheta_d=[0; 0]; % Desired velocity (Derivative of theta_d)
            ddtheta_d=[0; 0];
                
        end
        
        theta= x(1:2,1);
        dtheta= x(3:4,1);
        
        global Mmat Cmat
        Mmat =[I1+((l^2)*m2*(cos(x(2))^2))/4, 0;...
            0, I2+((l^2)*m2)/4];
        Cmat = [-(l^2*sin(2*x(2))*(m2)*x(4))/8, -(l^2*sin(2*x(2))*(m2)*x(3))/8 ;...
            (l^2*sin(2*x(2))*(m2)*(x(3)^2))/8, 0];
        J1 = [-l*sin(x(1))*cos(x(2)), -l*cos(x(1))*sin(x(2));...
            l*cos(x(1))*cos(x(2)), -l*sin(x(1))*sin(x(2));...
            0, l*cos(x(2))];
        invM = inv(Mmat);
        invMC = invM*Cmat;
        
        tau = PDControl(theta_d,dtheta_d,ddtheta_d,theta,dtheta, J1);
        
        torque =[torque, tau];
        dx=zeros(4,1);
        dx(1) = x(3); %dtheta1
        dx(2) = x(4); %dtheta2
        dx(3:4) = -invMC*x(3:4) +invM*tau; % because ddot theta = -M^{-1}(C \dot Theta) + M^{-1} tau
    end

    function dx = throwArmODE(t,x)
        
        theta_d = pi/2 ;
        dtheta_d = wf ; % Desired velocity (Derivative of theta_d)
        ddtheta_d = 0 ;
        
        theta = x(1);
        dtheta = x(2);
        
        global M C G
        M = I2 + (l^2*m2)/4 ;
        C = 0;
        G = 0.5*m2*g*l*cos(x(1)) ;
        J = [-l*sin(x(1)); 0; l*cos(x(1))];
        invM = 1/M;
        invMC = invM*C;
        
        if x(1)<= pi/2 
            tau = computeTorque(theta_d,dtheta_d,ddtheta_d,theta,dtheta, J);
%             omega=[omega;dtheta];
        elseif x(1)> pi/2
            tau = complianceControl(theta,dtheta, J);
            omega=[omega;dtheta];
        end
        
        
        tau_append = [0;tau] ;
        torque =[torque, tau_append];
        dx = zeros(2,1);
        
        dx(1)=x(2);
        dx(2)=-invMC + invM*tau;
        
        
    end

 function tau = PDControl(theta_d,dtheta_d,ddtheta_d,theta,dtheta, J)
     
        Kp=10*eye(2);
        Kv=10*eye(2);
        e=theta_d-theta; % position error
        de = dtheta_d - dtheta; % velocity error
        he =[0; 0; -0.1*9.81];
        tau = Kp*e + Kv*de -J'*he;
 end

function tau = computeTorque(theta_d, dtheta_d,ddtheta_d, theta, dtheta, J)
        global M C
        Kp = 25;
        Kv = 25;
        e = theta_d-theta ; % position error
        de = dtheta_d - dtheta; % velocity error
        he = [0;0;-0.1*9.81];
        tau= M*(Kp*e + Kv*de) + C*dtheta_d + M*ddtheta_d - J'*he;
end

function tau = complianceControl(theta, dtheta, J)
        Kpx = 20 ;
        Kdx = 40 ;
        Kp = Kpx*eye(3);
        Kv = Kdx*eye(3);
        e = pi -theta ; % position error
        de = 0 - dtheta; % velocity error
         xd = -l;
         xr = 0 ;
         fx = Kpx*(xd - xr) ;
        
         zd = 0;
         zr = l ;
         fz = Kdx*(zd - zr) ;
         he=[fx; 0; fz];
         tau= J'*Kp*J*e + J'*Kv*J*de - J'*he;
end

omega_actual=omega(1);
vel_actual=omega_actual*l;

figure
AxesH = axes('NextPlot', 'add');
time  = 0.46;
[x,y] = trajectory(vf, time);
plot(x, y, 'r', 'Parent', AxesH)
[x,y] = trajectory(vel_actual, time);
plot(x, y, 'b', 'Parent', AxesH)
xlabel('Distance in the XY Plane (m)');
ylabel('Z-Axis (m)');
legend('Desired Trajectory','Actual Trajectory') ;
box on;
grid on;
axis([0 1.5 0 1])


function [x,y]=trajectory(v,time)
x0=0;
y0=1;
g=9.81;
t=0:0.01:time;
x=x0+v*t;
y=y0-4.9*t.^2;
end

disp('Finish.');

end