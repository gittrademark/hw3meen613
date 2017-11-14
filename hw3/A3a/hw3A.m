% HW Assignment 3
% Sam Friedman, Benson Isaac, Mohamed Mohamed, Alexis Trevino
% 11/14/2017
% Script for dynamic analysis of System A 
clear
close all
%Initilize constants given in problem statement
L1 = 0.25;     % in m
L2 = 0.25;     % in m
M2 = 2;       % in kg
g  = 9.81;    % in m/(s^2)
w2 = M2*g;    % in N
phi0 = 45;
theta0 = 0;


%% Modal Analysis and linearization
%Create the matricies for linearized case
M = [  M2*L1^2      (M2*L1*L2)/2  ;...
      (M2*L1*L2)/2  (M2*L2^2)/3];
K = [  L1*w2         0            ;...
       0            (w2*L2)/2];

% Evaluate the eigenvalues, lambda, and corresponding eigenvectors, v.
[v,lambda] = eig(K,M);



%% Transient Analysis
%Choose time interval to evaluate dynamics and the time span.
tspan = 0:1/100:20;

%Choose initial condition for phi 
y0 = [phi0*pi/180 0 theta0*pi/180 0];
options = odeset('mass','M(t,y)');
[t,y]=ode113('indmot_ode',tspan,y0,options,L1,L2,M2,g,w2);

%Plot theta and phi over time
figure(1)
subplot(2,1,1)
plot(t,y(:,1)*180/pi);
grid on
hold on
xlabel('Time (s)');
ylabel('\phi (degrees)');
title('System A : Angle made by the string (\phi) for initial condition \phi_0 = 45^{\circ} with time')
ax = gca;
ax.XLim = [0 20];

subplot(2,1,2)
plot(t,y(:,3)*180/pi);
xlabel('Time (s)');
grid on
hold on
ylabel('\theta (degrees)');
title('System A : Angle made by the bar (\theta) for initial condition \phi_0 = 45^{\circ} with time')
ax = gca;
ax.XLim = [0 20];

% Calculation for tension
% Full, nonlinearized matrices based on MM.m and FF.m
X = zeros(size(y,1),2);
for i = 1:size(y,1)
    Mmatrix = [M2*L1^2                   , (M2*L2*L1/2)*cos(y(i,1)-y(i,3));
        (M2*L1*L2)/2*cos(y(i,3)-y(i,1)),  M2*L2^2/3                ];
    Fmatrix = [-w2*L1*sin(y(i,1))-(M2*L1*L2/2)*(y(i,4)^2)*sin(y(i,1)-y(i,3));
        (-w2*L2*sin(y(i,3)))/2+(M2*L1*L2/2)*(y(i,2)^2)*sin(y(i,1)-y(i,3))];
    X(i,:) = Mmatrix\Fmatrix;  % solve for phidt2, thetadt2
end
phi = y(:,1);
phidt = y(:,2);
phidt2 = X(:,1);
theta = y(:,3);
thetadt = y(:,4);
thetadt2 = X(:,2);

xgdt2 = -L1.*sin(phi).*(phidt.^2)+L1.*cos(phi).*phidt2 - ...
         0.5*L2.*sin(theta).*(thetadt.^2)+0.5*L2.*cos(theta).*thetadt2;
ygdt2 = L1.*cos(phi).*(phidt.^2)+L1.*sin(phi).*phidt2 + ...
        0.5*L2.*cos(theta).*(thetadt.^2)+0.5*L2.*sin(theta).*thetadt2;

tension = M2*sqrt(xgdt2.^2 + (ygdt2+g).^2);

figure(2)
plot(t,tension);
hold on
xlabel('Time (s)')
ylabel('Tension (N)')
suptitle('System A: Tension in the string when released at \phi_0 = 45^{\circ} ')

%% Frequency response
n = 1024;
y1 = fft(y(:,1),n);
m = abs(y1);
p = unwrap(angle(y1));
f = (0:length(y1)-1)*100/length(y1);

y2 = fft(y(:,3),n);
m2 = abs(y2);
p2 = unwrap(angle(y2));
f2 = (0:length(y2)-1)*100/length(y2);

figure(4)
subplot(2,1,1)
plot(f,m)
hold on
plot(f2,m2)
xlabel('Frequency (Hz)');
ylabel('Magnitude');
legend('\phi','\theta')
title('System A: Magnitude vs Frequency')
ax = gca;
ax.XLim = [0 20];


subplot(2,1,2)
plot(f,p)
hold on
plot(f2,p2)
xlabel('Frequency (Hz)');
ylabel('Phase(radians)');
legend('\phi','\theta')
title('System A: Phase vs Frequency')
ax = gca;
ax.XLim = [0 20];

