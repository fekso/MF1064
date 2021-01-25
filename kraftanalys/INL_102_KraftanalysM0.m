clc
clear all

%% Definiering av givna värden

L = 0.500;               %[m]
L1 = 0.657;              %[m]
L2 = 0.050;              %[m]
L3 = 0.670;              %[m]
L4x = 0.034;             %[m]
L4y = 0.044;             %[m]
h1 = 0.020;              %[m]
h2 = 0.050;              %[m]
h3 = 0.050;              %[m]
egenvikt = 'Vad väger en styrarm (0 eller 3)? '; % 3kg
m1 = input(egenvikt);
m2 = 10;                 %[kg]
m3 = 20;                 %[kg]
MC = 0;
g = 9.82;                %[m/s^2]
n = 2.5;                 %[varv per minut]
omega = 2*pi*n/60;       %[rad/s]
alpha  = [0:0.5:148.5]'*pi/180;
t = alpha/omega;

beta = 35*pi/180;     % Startgissning för vinkeln beta
gamma = 45*pi/180;    % Startgissning för vinkeln alfa
x = [beta; gamma];

%% Vi tar fram beta och gamma
for ii =1:length(alpha)
    
x = newtonrap(x,alpha(ii));     % Bestämning av beta & gamma med hjälp av Newton Raphson-funktionen

beta(ii) = x(1);
gamma(ii) = x(2);

end

beta = beta';
gamma = gamma';

 %% Analytisk och numerisk derivering av vinklar

for a = 1:length(alpha)
    xprim2 = analderiv(alpha(a),beta(a),gamma(a),L1,L2,L3,t,omega);
    betaprimanal(a) = xprim2(1);
    gammaprimanal(a) = xprim2(2);
end

xprim = numderiv(beta,gamma,alpha,t,omega);     % Matris med de numeriskt framtagna vinkelhastigheterna
xprim2 = [betaprimanal' gammaprimanal'];        % Matris med de analytiskt framtagna vinkelhastigheterna


%% Numeriska värden insatta i tvångsekvationerna

f1 = L1*cos(alpha)+L2*cos(gamma)-L3*cos(beta)-L4x;
f2 = L1*sin(alpha)-L2*sin(gamma)-L3*sin(beta)+L4y;

%% Masscentrum

mtot = 2*m1+m2+m3;

for i = 1:length(alpha);

mx1 = m1*L1*cos(alpha(i));
mx2 = m2*(L1*cos(alpha(i))+(L2-L)/2*cos(gamma(i)));
mx3 = m3*(L4x+L3/2*cos(beta(i)));
my1 = m1*L1*sin(alpha(i));
my2 = m2*(L1*sin(alpha(i))+(L-L2)/2*sin(gamma(i)));
my3 = m3*(L3/2*sin(beta(i))-L4y);

mxtot(i) = 1/(mtot)*(mx1+mx2+mx3);
mytot(i) = 1/(mtot)*(my1+my2+my3);

end

%% Plot av masscentrum

figure(1) 
plot(0,0,'o')
hold on
plot(mxtot,mytot)

xminT = -0.6;
xmaxT = 0.6;
yminT = -0.1;
ymaxT = 0.6;

set(findobj(gca,'type','line'), 'LineWidth', 3);
xlabel('x [m]')
ylabel('y [m]')
axis equal
xlim([xminT xmaxT]);
ylim([yminT ymaxT])
set(gca,'FontSize',16);
title('Masscentrums rörelse');
set(gca,'FontSize',36);
legend('Origo');
grid on

figure (2)

plot(alpha*180/pi,mxtot)
hold on
plot(alpha*180/pi,mytot)
set(findobj(gca,'type','line'), 'LineWidth', 3);
xlabel('\alpha[°]')
ylabel('[m]')
set(gca,'FontSize',16);
title('Masscentrums rörelse');
set(gca,'FontSize',36);
legend('x-koordinat','y-koordinat');
grid on

%% Jämviktsekvationer

    for j = 1:length(alpha)

    e1 = [1 0 -1 0 0 0 0 0 0];
    e2 = [0 1 0 1 0 0 0 0 0];
    e3 = [L2*sin(gamma(j)) L2*cos(gamma(j)) 0 0 0 0 0 0 0];
    e4 = [-1 0 0 0 1 0 0 0 0];
    e5 = [0 -1 0 0 0 1 0 0 0];
    e61 = [-L1*sin(alpha(j)) L1*cos(alpha(j)) 0 0 0 0 0 0 -1];;
    e7 = [0 0 1 0 0 0 1 0 0];
    e8 = [0 0 0 -1 0 0 0 1 0];
    e91 = [0 0 L3*sin(beta(j)) L3*cos(beta(j)) 0 0 0 0 0];

    AO = zeros(9,9);
    AO = [e1; e2; e3; e4; e5; e61; e7; e8; e91];
    
    b = [0 m2*g m2*g*(L2+L)/2*cos(gamma(j)) 0 2*m1*g -m1*g*L1*cos(alpha(j)) 0 m3*g -m3*g*L3/2*cos(beta(j))]';

    R1 = AO\b;

    rAX1(j) = R1(1);
    rAY1(j) = R1(2);
    rBX1(j) = R1(3);
    rBY1(j) = R1(4);
    rOX1(j) = R1(5);
    rOY1(j) = R1(6);
    rCX1(j) = R1(7);
    rCY1(j) = R1(8);
    MO(j) = R1(9);
    
    %% Validering 
    
    jmv1(j) = rAX1(j)-rBX1(j);
    jmv2(j) = rAY1(j)+rBY1(j)-m2*g;
    jmv3(j) = rAX1(j)*L2*sin(gamma(j))-rAY1(j)*L2*cos(gamma(j))-m2*g*(L1+L2)/2*cos(gamma(j));
    jmv4(j) = rOX1(j)-rAX1(j);
    jmv5(j) = rOY1(j)-rAY1(j)-2*m1*g;
    jmv6(j) = rAY1(j)*L1*cos(alpha(j))-rAX1(j)*L1*sin(alpha(j))+m1*g*L1*cos(alpha(j))-MO(j);
    jmv7(j) = rBX1(j)+rCX1(j);
    jmv8(j) = rCY1(j)-m3*g-rBY1(j);
    jmv9(j) = rBX1(j)*L3*sin(beta(j))+rBY1(j)*L3*cos(beta(j))+m3*g*L3/2*cos(beta(j))-MC;

    PO(j) = abs(MO(j)*omega); %effekt
    
 
    end
    
    %% Plot av fel
    
figure (3)
subplot(2,1,1);

plot(alpha*180/pi,jmv1)
xlabel('\alpha[°]')
ylabel('[Fel [N]')
title('Fel i jämviktsekvationer');
set(gca,'FontSize',12);
grid on
hold on
plot(alpha*180/pi,jmv2)
hold on
plot(alpha*180/pi,jmv4)
hold on
plot(alpha*180/pi,jmv5)
hold on
plot(alpha*180/pi,jmv7)
hold on
plot(alpha*180/pi,jmv8)
legend('Jämvikt 1','Jämvikt 2','Jämvikt 4','Jämvikt 5','Jämvikt 7','Jämvikt 8');
set(gca,'FontSize',24);
axis([0 150 -5*10^-14 5*10^-14])

subplot(2,1,2);
plot(alpha*180/pi,jmv3)
xlabel('\alpha[°]')
ylabel('Fel [Nm]')
title('Fel i jämviktsekvationer');
set(gca,'FontSize',12);
grid on
hold on
plot(alpha*180/pi,jmv6)
hold on
plot(alpha*180/pi,jmv9)
legend('Jämvikt 3','Jämvikt 6','Jämvikt 9');
set(gca,'FontSize',24);
axis([0 150 -10^-13 10^-13])
    
%% Resultanter

rAres = sqrt(rAX1.^2+rAY1.^2);
rBres = sqrt(rBX1.^2+rBY1.^2);
rOres = sqrt(rOX1.^2+rOY1.^2);
rCres = sqrt(rCX1.^2+rCY1.^2);


%% Plot över moment

figure (4)

plot(alpha*180/pi,MO)
set(findobj(gca,'type','line'), 'LineWidth', 3);
xlabel('\alpha[°]')
ylabel('[Nm]')
set(gca,'FontSize',16);
title('Moment i O av som funktion av \alpha');
set(gca,'FontSize',36);
grid on


%% Plot över effekt

figure (5)

plot(alpha*180/pi,PO)
set(findobj(gca,'type','line'), 'LineWidth', 3);
xlabel('\alpha[°]')
ylabel('[W]')
set(gca,'FontSize',16);
title('Drivande effekt');
set(gca,'FontSize',36);
grid on

%% Plot av reaktionskrafter

figure(6);

subplot(2,2,1); %rA
plot(alpha*180/pi,rAX1);
xlabel('\alpha [°]')
ylabel('Kraft [N]')
set(findobj(gca,'type','line'), 'LineWidth', 2);
set(gca,'FontSize',24);
hold on

plot(alpha*180/pi,rAY1);
set(findobj(gca,'type','line'), 'LineWidth', 2);
hold on

plot(alpha*180/pi,rAres);
set(findobj(gca,'type','line'), 'LineWidth', 2);

title({'\boldmath$r_A$'},'Interpreter','latex');
legend('x','y', 'resultant');
grid on
grid minor

subplot(2,2,2); %rB
plot(alpha*180/pi,rBX1);
xlabel('\alpha [°]')
ylabel('Kraft [N]')
set(findobj(gca,'type','line'), 'LineWidth', 2);
set(gca,'FontSize',24);
hold on

plot(alpha*180/pi,rBY1);
set(findobj(gca,'type','line'), 'LineWidth', 2);
hold on

plot(alpha*180/pi,rBres);
set(findobj(gca,'type','line'), 'LineWidth', 2);

title({'\boldmath$r_B$'},'Interpreter','latex');
legend('x','y', 'resultant');
grid on
grid minor

subplot(2,2,3); %rO
plot(alpha*180/pi,rOX1);
xlabel('\alpha [°]')
ylabel('Kraft [N]')
set(findobj(gca,'type','line'), 'LineWidth', 2);
set(gca,'FontSize',24);
hold on

plot(alpha*180/pi,rOY1);
set(findobj(gca,'type','line'), 'LineWidth', 2);
hold on

plot(alpha*180/pi,rOres);
set(findobj(gca,'type','line'), 'LineWidth', 2);

title({'\boldmath$r_O$'},'Interpreter','latex');
legend('x','y', 'resultant');
grid on
grid minor

subplot(2,2,4); %rC
plot(alpha*180/pi,rCX1);
xlabel('\alpha [°]')
ylabel('Kraft [N]')
set(findobj(gca,'type','line'), 'LineWidth', 2);
set(gca,'FontSize',24);
hold on

plot(alpha*180/pi,rCY1);
set(findobj(gca,'type','line'), 'LineWidth', 2);
hold on

plot(alpha*180/pi,rCres);
set(findobj(gca,'type','line'), 'LineWidth', 2);

title({'\boldmath$r_C$'},'Interpreter','latex');
legend('x','y', 'resultant');
grid on
grid minor


%% Energi

Ep = mtot*g*(mytot(end)-mytot(1));  %Joule/Newtonmeter
Emek = trapz(alpha,MO);             %Joule/Newtonmeter
Edelta = abs(Ep-Emek);              %Joule/Newtonmeter


%% Frågor

maxmoment = max(MO)
maxeffekt = max(PO)
Reaktionmax_O = max(rOres)
Reaktionmax_A = max(rAres)
Reaktionmax_B = max(rBres)
Reaktionmax_C = max(rCres)