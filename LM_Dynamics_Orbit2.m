%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Projeto Final do curso de Engenharia de Controle e Automação
% Universidade: CEFET - RJ/Uned NI
% Aluna: Laís Lima - Matrícula: 1620368ECAN
% Professor orientador: Mauro Vasconcellos
% Referência principal: Artigo "Three-Dimensional Trajectory Optimization of Soft Lunar Landings from the Parking Orbit with Considerations of the Landing Site" escrito por Bong-Gyun Park and Min-Jea Tahk (2011)
% Script: LM_Dynamics_Orbit2.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Este código apresenta a simulação da segunda órbita lunar considerada (aproximada por inspeção da referência principal citada anteriormente e testes prévios feitos com limitações de valores para os parâmetros r
% e vr nos momentos finais da missão de aterrissagem) para verificação da dinâmica do sistema do módulo lunar.

%% Seção 1: Análise em malha aberta
% Nesta seção do código, buscam-se pontos de operação para linearização do sistema 3D de um Módulo Lunar, através de aproximações das curvas obtidas no artigo de referência, cuja abordagem é um controle ótimo em malha aberta.

% Limpando o espaço de trabalho, a janela de comando e área de figuras do programa
clc; close all; clear;

% Parâmetros iniciais e finais 
Tf=3435.0; Tf2=3120.0; Tf1=3100.0; Tf0=2880.0; rr0=1844.2; %Rf=1737.4311589983485646371264010668;  Vvrf=-0.03897785431258966870604609766815; 

% Representações do sistema dinâmico do módulo lunar nas três fases do problema de aterrissagem
% De-Orbit Phase (Fase de queima de saída de órbita)
[t1,y1] = ode45(@(t,x)DeOrbitFunction(t,x,0,(-205+(25/7.37)*t)*(pi/180),1),[0 7.37],[rr0;-142.066*(pi/180);-20.24*(pi/180);0;1627.8e-3;-60.139e-3;600]);
                                                                                 
% Transfer Phase (Fase de transferência de órbita)
[t2,y2] = ode45(@(t,x)TransferOrbitFunction(t,x,0,-180*(pi/180),0),[7.37 7.37+Tf0],y1(end,:)');

% Powered Descent Phase (Fase de descida motorizada)
% Esta fase foi subdividida em três partes para melhor ajustar a aproximação das entradas (alpha, betha, k)
[t3,y3] = ode45(@PDFp1,[7.37+Tf0 Tf1],y2(end,:)'); % Powered Descent Phase - part 1                                                             
[t4,y4] = ode45(@PDFp2,[Tf1 Tf2],y3(end,:)'); % Powered Descent Phase - part 2                                                              
[t5,y5] = ode45(@PDFp3,[Tf2 Tf],y4(end,:)'); % Powered Descent Phase - part 3                                                             

% Reunindo as respostas do sistema
t=[t1;t2;t3;t4;t5];
y=[y1;y2;y3;y4;y5];

%% Seção 2: Gráficos dos estados do sistema na solução em malha aberta
% Neste seção os resultados obtidos como resposta do sistema dinâmico são apresentados graficamente.

% Distância radial
figure
plot(t1,y1(:,1),'b-')
hold on
plot(t2,y2(:,1),'b-')
hold on
plot(t3,y3(:,1),'b-')
hold on
plot(t4,y4(:,1),'b-')
hold on
plot(t5,y5(:,1),'b-')
title('Estado 1 - Distância radial')
ylabel('r (km)')
xlabel('t(s)')

% Longitude
figure
plot(t1,y1(:,2)*180/pi,'r-')
hold on
plot(t2,y2(:,2)*180/pi,'r-')
hold on
plot(t3,y3(:,2)*180/pi,'r-')
hold on
plot(t4,y4(:,2)*180/pi,'r-')
hold on
plot(t5,y5(:,2)*180/pi,'r-')
title('Estado 2 - Longitude')
ylabel('theta (°)')
xlabel('t(s)')

% Latitude
figure
plot(t1,y1(:,3)*180/pi,'g-')
hold on
plot(t2,y2(:,3)*180/pi,'g-')
hold on
plot(t3,y3(:,3)*180/pi,'g-')
hold on
plot(t4,y4(:,3)*180/pi,'g-')
hold on
plot(t5,y5(:,3)*180/pi,'g-')
title('Estado 3 - Latitude')
ylabel('phi (°)')
xlabel('t(s)')

% Velocidade radial
figure
plot(t1,y1(:,4)*1000,'m-')
hold on
plot(t2,y2(:,4)*1000,'m-')
hold on
plot(t3,y3(:,4)*1000,'m-')
hold on
plot(t4,y4(:,4)*1000,'m-')
hold on
plot(t5,y5(:,4)*1000,'m-')
title('Estado 4 - Veloc. Radial')
ylabel('vr (m/s)')
xlabel('t(s)')

% Velocidade longitudinal
figure
plot(t1,y1(:,5)*1000,'y-')
hold on
plot(t2,y2(:,5)*1000,'y-')
hold on
plot(t3,y3(:,5)*1000,'y-')
hold on
plot(t4,y4(:,5)*1000,'y-')
hold on
plot(t5,y5(:,5)*1000,'y-')
title('Estado 5 - Veloc. Longitudinal')
ylabel('vtheta (m/s)')
xlabel('t(s)')

% Velocidade latitudinal
figure
plot(t1,y1(:,6)*1000,'c-')
hold on
plot(t2,y2(:,6)*1000,'c-')
hold on
plot(t3,y3(:,6)*1000,'c-')
hold on
plot(t4,y4(:,6)*1000,'c-')
hold on
plot(t5,y5(:,6)*1000,'c-')
title('Estado 6 - Veloc. Latitudinal')
ylabel('vphi (m/s)')
xlabel('t(s)')

% Massa
figure
plot(t1,y1(:,7),'k-')
hold on
plot(t2,y2(:,7),'k-')
hold on
plot(t3,y3(:,7),'k-')
hold on
plot(t4,y4(:,7),'k-')
hold on
plot(t5,y5(:,7),'k-')
title('Estado 7 - Massa')
ylabel('m (kg)')
xlabel('t(s)')

% Gráfico de fase (Latitude x Longitude)
figure
plot(y1(:,2)*180/pi,y1(:,3)*180/pi,'m-')
hold on
plot(y2(:,2)*180/pi,y2(:,3)*180/pi,'m-')
hold on
plot(y3(:,2)*180/pi,y3(:,3)*180/pi,'m-')
hold on
plot(y4(:,2)*180/pi,y4(:,3)*180/pi,'m-')
hold on
plot(y5(:,2)*180/pi,y5(:,3)*180/pi,'m-')
title('Gráfico de Fase')
xlabel('Longitude (°)')
ylabel('Latitude (°)')

% Altitude x Longitude x Latitude
figure
plot3(y(:,3)*180/pi,y(:,2)*180/pi,y(:,1)-1737.4)
xlabel('Latitude (°)')
ylabel('Longitude (°)')
zlabel('Altitude (km)')
title('Gráfico Altitude x Longitude x Latitude')

%% Seção 3: Funções representativas da dinâmica do sistema 3D em cada fase do problema de aterrissagem
% Finalmente, esta seção contém todas as funções correspondentes ao comportamento do sistema em cada fase do problema de aterrissagem.

% De-Orbit Phase (Fase de de-órbita)
function dydt = DeOrbitFunction(t,x,alpha,betha,k)
g0=9.81e-3; % Gravidade da Terra em Km/s^2
Isp=316; % Impulso Específico em segundos
Tmax=1700e-3; % Força de Propulsão Máxima em Kg*Km/s^2
mi=4902.78; % Parâmetro Gravitacional Padrão em km^3/s^2
omega=2.6632e-6; % Velocidade Angular da Lua em rad/s.
dydt=[x(4); 
      x(5)/(x(1)*cos(x(3))); 
      x(6)/x(1); 
      ((-Tmax*k)/x(7))*sin(betha) - (mi)/(x(1)^2) + (x(6)^2)/x(1) + (x(5)^2)/x(1) + x(1)*(omega^2)*cos(x(3)) + 2*omega*x(5)*cos(x(3));
      ((Tmax*k)/x(7))*cos(betha)*cos(alpha) - (x(4)*x(5))/x(1) + (x(5)*x(6)*sin(x(3)))/(x(1)*cos(x(3))) + 2*omega*x(6)*sin(x(3)) - 2*omega*x(4)*cos(x(3));
      ((Tmax*k)/x(7))*cos(betha)*sin(alpha) - (x(4)*x(6))/x(1) - ((x(5)^2)*sin(x(3)))/(x(1)*cos(x(3))) - x(1)*(omega^2)*sin(x(3))*cos(x(3)) - 2*omega*x(5)*sin(x(3));
      -(Tmax*k)/(Isp*g0)];
end

% Transfer Orbit Phase (Fase de transferência de órbita)
function dydt = TransferOrbitFunction(t,x,alpha,betha,k)
g0=9.81e-3; % Gravidade da Terra em Km/s^2
Isp=316; % Impulso Específico em segundos
Tmax=1700e-3; % Força de Propulsão Máxima em Kg*Km/s^2
mi=4902.78; % Parâmetro Gravitacional Padrão em km^3/s^2
omega=2.6632e-6; % Velocidade Angular da Lua em rad/s.
dydt=[x(4); 
      x(5)/(x(1)*cos(x(3))); 
      x(6)/x(1); 
      ((-Tmax*k)/x(7))*sin(betha) - (mi)/(x(1)^2) + ((x(6)^2))/x(1) + ((x(5))^2)/x(1) + x(1)*(omega^2)*cos(x(3)) + 2*omega*x(5)*cos(x(3));
      ((Tmax*k)/x(7))*cos(betha)*cos(alpha) - (x(4)*x(5))/x(1) + (x(5)*x(6)*sin(x(3)))/(x(1)*cos(x(3))) + 2*omega*x(6)*sin(x(3)) - 2*omega*x(4)*cos(x(3));
      ((Tmax*k)/x(7))*cos(betha)*sin(alpha) - (x(4)*x(6))/x(1) - ((x(5)^2)*sin(x(3)))/(x(1)*cos(x(3))) - x(1)*(omega^2)*sin(x(3))*cos(x(3)) - 2*omega*x(5)*sin(x(3));
      -(Tmax*k)/(Isp*g0)];
end

%%% Powered Descent Phase (Fase de descida motorizada)
% Powered Descent Phase - part 1
function dydt = PDFp1(t,x)
betha=(-180+30*(t-2880.37)/(3510-2880.37))*pi/180;
alpha=(6.5*(t-2880.37)/219.63)*pi/180;
k=1;
g0=9.81e-3; % Gravidade da Terra em Km/s^2
Isp=316; % Impulso Específico em segundos
Tmax=1700e-3; % Força de Propulsão Máxima em Kg*Km/s^2
mi=4902.78; % Parâmetro Gravitacional Padrão em km^3/s^2
omega=2.6632e-6; % Velocidade Angular da Lua em rad/s.
dydt=[x(4); 
      x(5)/(x(1)*cos(x(3))); 
      x(6)/x(1); 
      ((-Tmax*k)/x(7))*sin(betha) - (mi)/(x(1)^2) + ((x(6)^2))/x(1) + ((x(5))^2)/x(1) + x(1)*(omega^2)*cos(x(3)) + 2*omega*x(5)*cos(x(3));
      ((Tmax*k)/x(7))*cos(betha)*cos(alpha) - (x(4)*x(5))/x(1) + (x(5)*x(6)*sin(x(3)))/(x(1)*cos(x(3))) + 2*omega*x(6)*sin(x(3)) - 2*omega*x(4)*cos(x(3));
      ((Tmax*k)/x(7))*cos(betha)*sin(alpha) - (x(4)*x(6))/x(1) - ((x(5)^2)*sin(x(3)))/(x(1)*cos(x(3))) - x(1)*(omega^2)*sin(x(3))*cos(x(3)) - 2*omega*x(5)*sin(x(3));
      -(Tmax*k)/(Isp*g0)];
end

% Powered Descent Phase - part 2
function dydt = PDFp2(t,x)
betha=(-180+30*(t-2880.37)/(3510-2880.37))*pi/180;
alpha=(6.5-0.5*(t-3100)/300)*pi/180;
k=1;
g0=9.81e-3; % Gravidade da Terra em Km/s^2
Isp=316; % Impulso Específico em segundos
Tmax=1700e-3; % Força de Propulsão Máxima em Kg*Km/s^2
mi=4902.78; % Parâmetro Gravitacional Padrão em km^3/s^2
omega=2.6632e-6; % Velocidade Angular da Lua em rad/s.
dydt=[x(4); 
      x(5)/(x(1)*cos(x(3))); 
      x(6)/x(1); 
      ((-Tmax*k)/x(7))*sin(betha) - (mi)/(x(1)^2) + ((x(6)^2))/x(1) + ((x(5))^2)/x(1) + x(1)*(omega^2)*cos(x(3)) + 2*omega*x(5)*cos(x(3));
      ((Tmax*k)/x(7))*cos(betha)*cos(alpha) - (x(4)*x(5))/x(1) + (x(5)*x(6)*sin(x(3)))/(x(1)*cos(x(3))) + 2*omega*x(6)*sin(x(3)) - 2*omega*x(4)*cos(x(3));
      ((Tmax*k)/x(7))*cos(betha)*sin(alpha) - (x(4)*x(6))/x(1) - ((x(5)^2)*sin(x(3)))/(x(1)*cos(x(3))) - x(1)*(omega^2)*sin(x(3))*cos(x(3)) - 2*omega*x(5)*sin(x(3));
      -(Tmax*k)/(Isp*g0)];
end

% Powered Descent Phase - part 3
function dydt = PDFp3(t,x)
betha=(-180+30*(t-2880.37)/(3510-2880.37))*pi/180;
alpha=(6+3*(t-3400)/110)*pi/180;
k=1;
g0=9.81e-3; % Gravidade da Terra em Km/s^2
Isp=316; % Impulso Específico em segundos
Tmax=1700e-3; % Força de Propulsão Máxima em Kg*Km/s^2
mi=4902.78; % Parâmetro Gravitacional Padrão em km^3/s^2
omega=2.6632e-6; % Velocidade Angular da Lua em rad/s.
dydt=[x(4); 
      x(5)/(x(1)*cos(x(3))); 
      x(6)/x(1); 
      ((-Tmax*k)/x(7))*sin(betha) - (mi)/(x(1)^2) + ((x(6)^2))/x(1) + ((x(5))^2)/x(1) + x(1)*(omega^2)*cos(x(3)) + 2*omega*x(5)*cos(x(3));
      ((Tmax*k)/x(7))*cos(betha)*cos(alpha) - (x(4)*x(5))/x(1) + (x(5)*x(6)*sin(x(3)))/(x(1)*cos(x(3))) + 2*omega*x(6)*sin(x(3)) - 2*omega*x(4)*cos(x(3));
      ((Tmax*k)/x(7))*cos(betha)*sin(alpha) - (x(4)*x(6))/x(1) - ((x(5)^2)*sin(x(3)))/(x(1)*cos(x(3))) - x(1)*(omega^2)*sin(x(3))*cos(x(3)) - 2*omega*x(5)*sin(x(3));
      -(Tmax*k)/(Isp*g0)];
end

