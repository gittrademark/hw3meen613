% HW Assignment 3
% Sam Friedman, Benson Isaac, Mohamed Mohamed, Alexis Trevino
% 11/14/2017
% Script for dynamic analysis of System B
l=0.25; L=0.25; m=2; g=9.81; w=m*g; d=0.3;
M2=m;
L1=l;
L2=L;
w2=w;
tspan= [0:1/100:20];
options=odeset('mass','M(t,y)');

% Initial condition for phi and corresponding theta
phi_init = 45*pi/180;
theta_init = acos((d-l*cos(phi_init))/L);

% Transient analysis
y0=[phi_init;0;theta_init;0;0];
[t,y]=ode113('indmot_ode1',tspan,y0,options,m,l,L);

% Plot for theta and phi
figure(1)
subplot(2,1,1)
plot(t,y(:,1)*180/pi);
grid on
xlabel('Time (s)');
ylabel('\phi (degrees)');
title('System B : Angle made by the string (\phi) for initial condition \phi_0 = 45^{\circ} with time')

subplot(2,1,2)
plot(t,y(:,3)*180/pi);
grid on
xlabel('Time (s)');
ylabel('\theta (degrees)');
title('System B : Angle made by the bar (\theta) for initial condition \phi_0 = 45^{\circ} with time')

% Plot for tension and normal force with time
% Full, nonlinearized matrices based on MM.m and FF.m
X = zeros(size(y,1),3);
for i = 1:size(y,1)
    Mmatrix = [M2*L1^2                   , (M2*L2*L1/2)*cos(y(i,1)-y(i,3)), L1*sin(y(i,1));
        (M2*L1*L2)/2*cos(y(i,3)-y(i,1)),  M2*L2^2/3, L2*sin(y(i,3))               ;
        L1*sin(y(i,1)), L2*sin(y(i,3)), 0];
    Fmatrix = [-w2*L1*sin(y(i,1))-(M2*L1*L2/2)*(y(i,4)^2)*sin(y(i,1)-y(i,3));
        (-w2*L2*sin(y(i,3)))/2+(M2*L1*L2/2)*(y(i,2)^2)*sin(y(i,1)-y(i,3));
        -L1*y(i,2)^2*cos(y(i,1))-L2*y(i,4)^2*cos(y(i,3)) ];
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

tension = -M2.*xgdt2./sin(phi); 
normalf = M2.*(ygdt2+g) - tension.*cos(phi); %also normalf = -X(:,3)
figure(5)
plot(t,tension);
hold on
plot(t,normalf)
xlabel('Time(s)')
ylabel('Force(N)')
legend('Tension','Normal')
suptitle('System B: Tension and Normal forces in the string when released at \phi_0 = 45^{\circ} ')

%% Frequency response
% Frequency response for phi
n = 512;
y1 = fft(y(:,1),n);
m = abs(y1);
m(1) =0;
p = unwrap(angle(y1));
f = (0:length(y1)-1)*100/length(y1);

% Frequency response for theta
y2 = fft(y(:,3),n);
m2 = abs(y2);
m2(1) =0;
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
title('System B: Magnitude vs Frequency')
ax = gca;
ax.XLim = [0 10];


subplot(2,1,2)
plot(f,p)
hold on
plot(f2,p2)
xlabel('Frequency (Hz)');
ylabel('Phase(radians)');
legend('\phi','\theta')
title('System B: Phase vs Frequency')
ax = gca;
ax.XLim = [0 10];