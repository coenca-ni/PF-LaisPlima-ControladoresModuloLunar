%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Projeto Final do curso de Engenharia de Controle e Automação
% Universidade: CEFET - RJ/Uned NI
% Aluna: Laís Lima - Matrícula: 1620368ECAN
% Professor orientador: Mauro Vasconcellos
% Referência principal: Artigo "Three-Dimensional Trajectory Optimization of Soft Lunar Landings from the Parking Orbit with Considerations of the Landing Site" escrito por Bong-Gyun Park and Min-Jea Tahk (2011)
% Script: LM_PureDynamics.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Este código deve ser rodado em primeira instância.

%% Seção 1: Análise em malha aberta
% Nesta seção do código, buscam-se pontos de operação para linearização do sistema 3D de um Módulo Lunar, através de aproximações das curvas obtidas no artigo de referência, cuja abordagem é um controle ótimo em malha aberta.

% Limpando o espaço de trabalho, a janela de comando e área de figuras do programa
clc; close all; clear all;

% Parâmetros iniciais e finais para a simulação do sistema dinâmico (ÓTIMOS)
Tff=3400; Tf=3350; Tf2=3340; Tf1=2926; Tf0=2862; rr0=1837.4; %Rf=1737.3979609789;  Vvrf=0,844185604183028 m/s; 

%Intervalos temporais referentes a cada fase - Abordagem com o comando linspace
span1=linspace(0,7.37,7.37); %old: [0 7.37]
span2=linspace(7.37, 7.37+Tf0,Tf0); %old: [7.37 7.37+Tf0]
spanp1=linspace(7.37+Tf0,Tf1,Tf1-Tf0-7.37); %old: [7.37+Tf0 Tf1]
spanp2=linspace(Tf1,Tf2,Tf2-Tf1); %old: [Tf1 Tf2]
spanp3=linspace(Tf2,Tf,Tf-Tf2); %old: [Tf2 Tf]
spanp4=linspace(Tf,Tff,Tff-Tf); %old: [Tf Tff]

%Intervalos temporais referentes a cada fase - Abordagem antiga
% span1=[0 7.37];
% span2=[7.37 7.37+Tf0];
% spanp1=[7.37+Tf0 Tf1];
% spanp2=[Tf1 Tf2];
% spanp3=[Tf2 Tf];
% spanp4=[Tf Tff];

% Representações do sistema dinâmico do módulo lunar nas três fases do problema de aterrissagem
% De-Orbit Phase (Fase de queima de saída de órbita)
[t1,y1] = ode45(@(t,x)DeOrbitFunction(t,x,0,(-205+(25/7.37)*t)*(pi/180),1),span1,[rr0;-142.066*(pi/180);-20.24*(pi/180);0;1627.8e-3;-60.139e-3;600]);
                                                                                 
% Transfer Phase (Fase de transferência de órbita)
[t2,y2] = ode45(@(t,x)TransferOrbitFunction(t,x,0,-180*(pi/180),0),span2,y1(end,:)');

% Powered Descent Phase (Fase de descida motorizada) - Terceira abordagem
% Esta fase foi subdividida em quatro partes para melhor ajustar a aproximação das entradas (alpha, betha, k)
[t3,y3] = ode45(@(t,x)PDFp1(t,x,Tf0,Tf1,Tf),spanp1,y2(end,:)'); % Powered Descent Phase - part 1                                                             
[t4,y4] = ode45(@(t,x)PDFp2(t,x,Tf0,Tf1,Tf2,Tf),spanp2,y3(end,:)'); % Powered Descent Phase - part 2                                                              
[t5,y5] = ode45(@(t,x)PDFp3(t,x,Tf0,Tf2,Tf),spanp3,y4(end,:)'); % Powered Descent Phase - part 3                                                             
[t6,y6] = ode45(@(t,x)PDFp4(t,x),spanp4,y5(end,:)');

% Reunindo as respostas do sistema
t=[t1;t2;t3;t4;t5;t6];
y=[y1;y2;y3;y4;y5;y6];

%% Seção 2: Truncamento da matriz y para plotagem dos gráficos dos estados

% Código para truncar a matriz y no instante em que r=1737.3979:
linhadeparada=0;
tempodeparada=0;
rfinal=0;
for i=1:1:length(t)
    if(t(i)==3331.980629539951678452780470252)
    linhadeparada=i;
    tempodeparada=vpa(t(i));
    rfinal=vpa(y(i,4));
    end
end

Y_PLOT=[];
T_PLOT=[];
contador=0;

for i=1:1:length(t)
    if(t(i)<=tempodeparada)
        T_PLOT=[T_PLOT; t(i);];
        Y_PLOT=[Y_PLOT; y(i,:);];
    end
end

% De acordo com o que já observamos, o sistema é controlável nas primeira e
% última fases. A última fase para nós compreende os resultados do instante Tf0 ao
% tempodeparada! Para saber quantos pontos estavam compreendidos neste trecho, fiz o seguinte laço e conclui que são 468 pontos: 
% t_testando=[];
% y_testando=[];
% for i=1:1:length(t)
%     if(t(i)>=Tf0 && t(i)<=tempodeparada)
%         t_testando=[t_testando; t(i);];
%         y_testando=[y_testando; y(i,:);];
%     end
% end

%% Seção 3: Gráficos dos estados "truncados", desconsiderando visualmente a fase adicional PDFp4

figure
plot(T_PLOT, Y_PLOT(:,1))
legend('Distância Radial (km)')
figure
plot(T_PLOT, Y_PLOT(:,2)*180/pi)
legend('Longitude (°)')
figure
plot(T_PLOT, Y_PLOT(:,3)*180/pi)
legend('Latitude (°)')
figure
plot(T_PLOT, Y_PLOT(:,4)*1000)
legend('Velocidade Radial (m/s)')
figure
plot(T_PLOT, Y_PLOT(:,5)*1000)
legend('Velocidade Longitudinal (m/s)')
figure
plot(T_PLOT, Y_PLOT(:,6)*1000)
legend('Velocidade Latitudinal (m/s)')
figure
plot(T_PLOT, Y_PLOT(:,7))
legend('Massa (kg)')


%% Seção 4: Funções representativas da dinâmica do sistema 3D em cada fase do problema de aterrissagem
% De-Orbit Phase
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

% Transfer Orbit Phase
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

%%% Powered Descent Phase - Terceira abordagem
% Powered Descent Phase - part 1
function dydt = PDFp1(t,x,Tf0,Tf1,Tf)
betha=(-180+30*(t-Tf0)/(Tf-Tf0))*pi/180;
alpha=(6.5*(t-Tf0)/(Tf1-Tf0))*pi/180;
k=1;
%k=0.0025*t-7.0093; %Aproximei pegando 2 pontos do gráfico (ver bloco de notas das condições iniciais)
g0=9.81e-3; % gravidade da Terra em Km/s^2
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
function dydt = PDFp2(t,x,Tf0,Tf1,Tf2,Tf)
betha=(-180+30*(t-Tf0)/(Tf-Tf0))*pi/180;
alpha=(6.5-0.5*(t-Tf1)/(Tf2-Tf1))*pi/180;
k=1;
g0=9.81e-3; % gravidade da Terra em Km/s^2
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
function dydt = PDFp3(t,x,Tf0,Tf2,Tf)
betha=(-180+30*(t-Tf0)/(Tf-Tf0))*pi/180;
alpha=(6+3*(t-Tf2)/(Tf-Tf2))*pi/180;
k=1;
g0=9.81e-3; % km/s^2
Isp=316; % segundos
Tmax=1700e-3; % kg*km/s^2
mi=4902.78; % km^3/s^2
omega=2.6632e-6; % rad/s.
dydt=[x(4); 
      x(5)/(x(1)*cos(x(3))); 
      x(6)/x(1); 
      ((-Tmax*k)/x(7))*sin(betha) - (mi)/(x(1)^2) + ((x(6)^2))/x(1) + ((x(5))^2)/x(1) + x(1)*(omega^2)*cos(x(3)) + 2*omega*x(5)*cos(x(3));
      ((Tmax*k)/x(7))*cos(betha)*cos(alpha) - (x(4)*x(5))/x(1) + (x(5)*x(6)*sin(x(3)))/(x(1)*cos(x(3))) + 2*omega*x(6)*sin(x(3)) - 2*omega*x(4)*cos(x(3));
      ((Tmax*k)/x(7))*cos(betha)*sin(alpha) - (x(4)*x(6))/x(1) - ((x(5)^2)*sin(x(3)))/(x(1)*cos(x(3))) - x(1)*(omega^2)*sin(x(3))*cos(x(3)) - 2*omega*x(5)*sin(x(3));
      -(Tmax*k)/(Isp*g0)];
end

% Powered Descent Phase - part 4
function dydt = PDFp4(t,x)
betha=0;
alpha=0;
k=0;
g0=9.81e-3; % km/s^2
Isp=316; % segundos
Tmax=1700e-3; % kg*km/s^2
mi=4902.78; % km^3/s^2
omega=2.6632e-6; % rad/s.
dydt=[x(4); 
      x(5)/(x(1)*cos(x(3))); 
      x(6)/x(1); 
      ((-Tmax*k)/x(7))*sin(betha) - (mi)/(x(1)^2) + ((x(6)^2))/x(1) + ((x(5))^2)/x(1) + x(1)*(omega^2)*cos(x(3)) + 2*omega*x(5)*cos(x(3));
      ((Tmax*k)/x(7))*cos(betha)*cos(alpha) - (x(4)*x(5))/x(1) + (x(5)*x(6)*sin(x(3)))/(x(1)*cos(x(3))) + 2*omega*x(6)*sin(x(3)) - 2*omega*x(4)*cos(x(3));
      ((Tmax*k)/x(7))*cos(betha)*sin(alpha) - (x(4)*x(6))/x(1) - ((x(5)^2)*sin(x(3)))/(x(1)*cos(x(3))) - x(1)*(omega^2)*sin(x(3))*cos(x(3)) - 2*omega*x(5)*sin(x(3));
      -(Tmax*k)/(Isp*g0)];
end
