clc
clear all

%% Definiering av givna värden

L = 600;              %[mm]
L1 = 657;             %[mm]
L2 = 50;             %[mm]
L3 = 670;             %[mm]
L4x = 35;            %[mm]
L4y = 44;             %[mm]
h1 = 20;              %[mm]
h2 = 50;              %[mm]
h3 = 50;              %[mm]
n = 2.5;              %[varv per minut]
omega = 2*pi*n/60;    %[rad/s]
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

%% Plot över vinkeln beta & gamma som funktion av alfa

figure(1)
plot(alpha*180/pi, beta*180/pi);
hold on
plot(alpha*180/pi, gamma*180/pi);
set(findobj(gca,'type','line'), 'LineWidth', 2)
set(gca,'FontSize',36);
xlabel('\alpha [°]')
title('Vinkeln \beta & \gamma som funktion av vinkeln \alpha');
legend('\beta','\gamma');
grid on
   
%% Numeriska värden insatta i tvångsekvationerna

f1 = L1*cos(alpha)+L2*cos(gamma)-L3*cos(beta)-L4x;
f2 = L1*sin(alpha)-L2*sin(gamma)-L3*sin(beta)+L4y;

figure(2)
plot(alpha*180/pi,f1)
hold on
plot(alpha*180/pi,f2)
set(findobj(gca,'type','line'), 'LineWidth', 2)
set(gca,'FontSize',36);
xlabel('\alpha [°]')
title('Tvångsekvationernas numeriska värden');
legend('\beta','\gamma');
grid on

%% Definiering av ortsvektorer

x0 = zeros(size(alpha));
xA = x0 + L1*cos(alpha);
y0 = zeros(size(alpha));
yA = y0 + L1*sin(alpha);

xB = zeros(size(alpha));
xB = xA + L2*cos(gamma);
yB = zeros(size(alpha));
yB = yA - L2*sin(gamma);

xP = zeros(size(alpha));
xP = xA - L*cos(gamma);
yP = zeros(size(alpha));
yP = yA + L*sin(gamma);

xC = zeros(size(alpha));
xC = xA + xB - L3*cos(beta);
yC = zeros(size(alpha));
yC = yA + yB - L3*sin(beta);

xT = zeros(size(alpha));
xT = xP + h2/2*sin(gamma);
yT = zeros(size(alpha));
yT = yP + h2/2*cos(gamma);

xmin = -1300;
ymin = -200;
xmax = 1000;
ymax = 1000;

%% Animering av konstruktionen och punkten T:s rörelse

figure(3);
for j = 1:length(alpha)
   figure(3)
   plot([x0(j),xA(j)],[y0(j),yA(j)]);
   hold on
   plot([xA(j),xP(j)],[yA(j),yP(j)]);
   hold on
   plot([xA(j),xB(j)],[yA(j),yB(j)]);
   hold on
   plot([xP(j),xT(j)],[yP(j),yT(j)]);
   hold on
   plot([xB(j),L4x],[yB(j),-L4y]);
   set(findobj(gca,'type','line'), 'LineWidth', 3);
   xlabel('x [mm]')
   ylabel('y [mm]')
   axis equal
   xlim([xmin xmax]);
   ylim([ymin xmax])
   hold off
   pause(0.01)
end

 %% Analytisk och numerisk derivering av vinklar

for a = 1:length(alpha)
    xprim2 = analderiv(alpha(a),beta(a),gamma(a),L1,L2,L3,t,omega);
    betaprimanal(a) = xprim2(1);
    gammaprimanal(a) = xprim2(2);
end

xprim = numderiv(beta,gamma,alpha,t,omega);     % Matris med de numeriskt framtagna vinkelhastigheterna
xprim2 = [betaprimanal' gammaprimanal'];        % Matris med de analytiskt framtagna vinkelhastigheterna


%% Skillnad mellan numerisk och analytisk derivering

figure (4);

subplot(1,2,1)
plot(alpha*180/pi,xprim(:,1));
hold on
plot(alpha*180/pi,xprim(:,2));
set(gca,'FontSize',16);
set(findobj(gca,'type','line'), 'LineWidth', 2)
title('Vinkelhastigheter för vinklarna \beta & \gamma bestämda numeriskt');
legend({'$\dot{\beta}(\alpha_i)-numerisk$', '$\dot{\gamma}(\alpha_i)-numerisk$'},'Interpreter','latex');
grid on

subplot(1,2,2)
plot(alpha*180/pi,xprim2(:,1));
hold on
plot(alpha*180/pi,xprim2(:,2));

set(gca,'FontSize',16);
set(findobj(gca,'type','line'), 'LineWidth', 2)
title('Vinkelhastigheter för vinklarna \beta & \gamma bestämda analytiskt');
legend({'$\dot{\beta}(\alpha_i)-analytsik$', '$\dot{\gamma}(\alpha_i)-analytisk$'},'Interpreter','latex');
grid on

%% Punkten T:s läge

figure(5) 
plot(0,0,'o')
hold on
plot(xT/1000,yT/1000)

xminT = -1.2;
xmaxT = 0.4;
yminT = -0.1;
ymaxT = 1;

set(findobj(gca,'type','line'), 'LineWidth', 3);
xlabel('x [m]')
ylabel('y [m]')
axis equal
xlim([xminT xmaxT]);
ylim([yminT ymaxT])
set(gca,'FontSize',16);
title('Punktens T:s rörelse');
set(gca,'FontSize',36);
legend('Origo');
grid on
   
figure(6)
plot(alpha*180/pi,yT/1000)
set(findobj(gca,'type','line'), 'LineWidth', 3);
xlabel('\alpha [°]')
ylabel('y [m]')
set(gca,'FontSize',36);
title('Punkten T:s höjd som funktion av \alpha');
grid on

%% Hastigheten för T analytiskt

for a = 1:length(alpha)
    dxTA(a) = -L1*omega*sin(alpha(a))+L*gammaprimanal(a)*sin(gamma(a))+h2/2*gammaprimanal(a)*cos(gamma(a));
    dyTA(a) = L1*omega*cos(alpha(a))+L*gammaprimanal(a)*cos(gamma(a))-h2/2*gammaprimanal(a)*sin(gamma(a));
    vTA(a) = sqrt(dxTA(a).^2+dyTA(a).^2);           % Hastighetsresultanten av x- och y-led
end

%% Hastigheten för T numeriskt

N = length(alpha);

for b = 2:(N-1)
    dt = (alpha(b+1)-alpha(b-1))/omega;
    dxTn = xT(b+1)-xT(b-1);
    dyTn = yT(b+1)-yT(b-1);
    dxTnum(b) = dxTn/dt;
    dyTnum(b) = dyTn/dt;
end

dxTnum(1) = NaN;
dyTnum(1) = NaN;
dxTnum(N) = NaN;
dyTnum(N) = NaN;

vTnum = sqrt(dxTnum.^2 + dyTnum.^2);                % Hastighetsresultanten av x- och y-led

%% Analytiska hastigheter

figure(7)

subplot(1,2,1);
plot(alpha*180/pi,dxTA/1000);
xlabel('\alpha [°]')
ylabel('v [m/s]')
set(findobj(gca,'type','line'), 'LineWidth', 2);
set(gca,'FontSize',24);
hold on

plot(alpha*180/pi,dyTA/1000);
set(findobj(gca,'type','line'), 'LineWidth', 2);
hold on

plot(alpha*180/pi,vTA/1000);
set(findobj(gca,'type','line'), 'LineWidth', 2);

title('Hastigheterna analytiskt bestämda');
legend('x','y', 'resultant');

%% Numeriska hastigheter

figure(7);

subplot(1,2,2);
plot(alpha*180/pi,dxTnum/1000);
xlabel('\alpha [°]')
ylabel('v [m/s]')
set(findobj(gca,'type','line'), 'LineWidth', 2);
set(gca,'FontSize',24);
hold on

plot(alpha*180/pi,dyTnum/1000);
set(findobj(gca,'type','line'), 'LineWidth', 2);
hold on

plot(alpha*180/pi,vTnum/1000);
set(findobj(gca,'type','line'), 'LineWidth', 2);

title('Hastigheterna numeriskt bestämda');
legend('x','y', 'resultant');
grid on

%%  Svar på frågorna

alphamax = alpha(end)*180/pi;
tmax = t(end);
betamax = max(beta)*180/pi;
Tmax = max(yT)/1000;
vTslutanal = vTA(end)/1000;
vTslutnum = vTnum(end-1)/1000;

Max_styrarm = num2str(alphamax)                     % Styrarmens största vinkel
Tid = num2str(tmax)                                 % Tiden det tar flr taket att fälla upp
Max_bakruta = num2str(betamax)                      % Största vinkeln bakrutan
Max_height_T = num2str(Tmax)                        % Högsta punkten taket antar under rörelsen
Hastighet_T_slut_analytisk = num2str(vTslutanal)    % Punkten T:s hastighet när taket är i uppfällt lägebestämt analytiskt 
Hastighet_T_slut_numerisk = num2str(vTslutnum)      % Punkten T:s hastighet när taket är i uppfällt lägebestämt numeriskt

