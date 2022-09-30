%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Projeto Final do curso de Engenharia de Controle e Automação
% % Universidade: CEFET - RJ/Uned NI
% % Aluna: Laís Lima - Matrícula: 1620368ECAN
% % Professor orientador: Mauro Vasconcellos
% % Referência principal: Artigo "Three-Dimensional Trajectory Optimization of Soft Lunar Landings from the Parking Orbit with Considerations of the Landing Site" escrito por Bong-Gyun Park and Min-Jea Tahk (2011)
% % Script: LQG.m
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Projetando o controlador LQG para o sistema 3D de um módulo lunar 

% Montando a matriz de saídas do sistema
% % C=eye(7);
%% C=[1 1e-5 1e-5 1e-5 1e-5 1e-5 1e-5;1e-5 1e-5 1e-5 1 1e-5 1e-5 1e-5];
% % C=[1 1 1 1 1 1 1];
% % C3=[1 1e-2 1e-2 1e-2 1e-2 1e-2 1e-2];
% C3=[1e3 0 0 0 0 0 0];
   close all;
   C2=[1 0 0 0 0 0 0; 
       0 0 0 0 0 0 0;   
       0 0 0 0 0 0 0]; %matriz de ajuste de dimensão
   C=[1e3 0 0 0 0 0 0; 
      0 0 0 0 0 0 0;   
      0 0 0 0 0 0 0;
      0 0 0 0 0 0 0; 
      0 0 0 0 0 0 0;
      0 0 0 0 0 0 0;
      0 0 0 0 0 0 0     
     ]; %considerando apenas o sensor de posição (r) 
                   
  dimC=size(C);
  linC=dimC(1,1);
  colC=dimC(1,2);

  D=0; 

% Definindo pesos componentes da matriz Q e R por tentativa e erro
%Q1=1; Q2=1; Q3=1; Q4=1; Q5=1; Q6=1; Q7=1; %nope
%Q1=.03; Q2=1; Q3=1; Q4=4.500; Q5=10; Q6=10; Q7=20; %nope
%Q1=1; Q2=1; Q3=1; Q4=10; Q5=10; Q6=10; Q7=.5; % começou a aproximar
%Q1=.03; Q2=1; Q3=1; Q4=4500; Q5=10; Q6=10; Q7=.5;
%Q1=.09; Q2=1; Q3=13.5; Q4=5; Q5=10; Q6=10; Q7=.5; % muito muito bom (mas vr
%ainda tá ruim no final)
% % % % % % %  Q1=.0485; Q1=.049;  Q1=.050515; Q1=0.00575; Q2=1; Q3=13.5; Q4=5e-3; Q5=10; Q6=10; Q7=.5; %último que foi para o texto
Q1=0.04942775; Q2=1; Q3=55; Q4=5e-3; Q5=10; Q6=10; Q7=.5;

%Q1=1000; Q2=1; Q3=1; Q4=1000; Q5=1; Q6=1; Q7=1;
Q=[Q1 0 0 0 0 0 0;
   0 Q2 0 0 0 0 0;
   0 0 Q3 0 0 0 0;
   0 0 0 Q4 0 0 0;
   0 0 0 0 Q5 0 0;
   0 0 0 0 0 Q6 0;
   0 0 0 0 0 0 Q7];
 Ru= diag([1,1,500]);
%Ru= diag([1,1,1]);

% Informações sobre valores válidos obtidos do código LM_DynamicsOptimized %%%%%
x_valido=[R_valido THETA_valido PHI_valido VR_valido VTHETA_valido VPHI_valido M_valido];

%%% Covariâncias dos Ruídos do Processo e de Medição %%%%%%%
W=.0100*eye(7,7); % Covariância da perturbação - Qk
V=.0100*eye(linC); % Covariância do ruído - Rk
Qk=W;
Rk=V;

%% Seção 2: Calculando as leis de controle 

QI=.100*eye(linC); % # linC saídas
QXU=blkdiag(Q,Ru); % matrizes de ponderação
QWV=blkdiag(Qk,Rk); % matrizes de covariância dos ruídos


inc=1; % inc=10
contador_t=-inc+1;
ii=-7*inc+1; 
REF=[]; TREF=[]; Aff=[];

R=[]; THETA=[]; PHI=[]; VR=[]; RKVK=[];
VTHETA=[]; VPHI=[]; M=[]; Xe=[]; W=[];

%%% Parâmetros e Condições Iniciais %%%
m0=600; % massa em Kg
mk=m0;
g0=9.81; % gravidade da Terra em m/s
Isp=316; % Impulso Específico em segundos
Tmax=1700; % Força de Propulsão Máxima em Newtons
r0=1837.4e3; % Posição Radial Inicial em m
um=4902.78e9; % Parâmetro Gravitacional Padrão em m^3/s^2
w=2.6632e-6; % Velocidade Angular da Lua em rad/s.
rk=r0;
theta0=-142.066; % Posição Longitudinal Inicial em graus
thetak=theta0*pi/180; % em radianos 
phi0=-20.24; % Posição Latitudinal Inicial em graus
phik=phi0*pi/180; % em radianos
vr0=0; % Velocidade Radial Inicial em m/s
vrk=vr0;
vtheta0=1627.8; % Velocidade Longitudinal Inicial em m/s
vthetak=vtheta0;
vphi0=-60.139; % Velocidade Latitudinal Inicial em m/s
vphik=vphi0;
xki=zeros(linC);
%%% estados estimados iniciais %%%
%xe0=[r0, ..., m0]; 
%xek=[rk thetak phik vrk vthetak vphik mk]'; % Estado Estimado
rke=rk;
thetake=thetak;
phike=phik;
vrke=vrk;
vthetake=vthetak;
vphike=vphik;
mke=mk;
%%% entradas iniciais %%%
alphak=0;
betak=-205*pi/180;
Kk=1;

% wxke=[rke thetake phike vrke vthetake vphike mke]';
% wxk=[rk thetak phik vrk vthetak vphik mk]';

% Condições Iniciais do Filtro de Kalman %%%%%%%%%%%%%
% xhatk-1=x(0) condição inicial
 Pk=Qk;

%%%%%%%%%%%%% Discretização de Euler %%%%%%%%%%%%%%%%

%for j=1:round((length(A_valida)/(7*inc)))-inc
for j=1:round(length(T_valido))-inc
    Ts=T_valido(j+1)-T_valido (j);  % intervalo de amostragem
    contador_t=contador_t+inc; 
    ii=ii+7*inc;
    KLQR=lqr(A_valida(ii:ii+6,:), B_valida(ii:ii+6,:),Q,Ru); 
%     SYS=ss(A_valida(ii:ii+6,:), B_valida(ii:ii+6,:),C,D);
%     [KLQG,INFO]=lqg(SYS,QXU,QWV,QI); 
%     
%     Af=KLQG.a;
%     Bf=KLQG.b;
%     Cf=KLQG.c;
%     Df=KLQG.d;

    Ak=A_valida(ii:ii+6,:)*Ts+eye(7);
    Bk=B_valida(ii:ii+6,:)*Ts;

    ref=C*x_valido(contador_t,:)';
    tref=T_valido(contador_t,:);
    REF=[REF ref(:,1)];
    TREF=[TREF tref]; T=TREF;

    % Equação de Estados da Planta %%%
    f1k=vrk;
    f2k=vthetak/(rk*cos(phik));
    f3k=vphik/rk;
    f4k=((-Tmax*Kk/mk)*sin(betak))+(-um/(rk^2))+(vphik^2/rk)+(vthetak^2/rk)+((rk*w^2)*(cos(phik)))+(2*w*vthetak*cos(phik));
    f5k=(((Tmax*Kk)/mk)*cos(betak)*cos(alphak))-(vthetak*vrk/rk)+((vphik*vthetak/rk)*sin(phik)/cos(phik))+(2*w*vphik*sin(phik))+(-2*w*vrk*cos(phik));
    f6k=(((Tmax*Kk)/mk)*cos(betak)*sin(alphak))+(-vphik*vrk/rk)+(((-vthetak^2)*sin(phik))/(cos(phik)*rk))+((-rk*w^2)*sin(phik)*cos(phik))+(-2*w*vthetak*sin(phik));
    f7k=(-Tmax*Kk/(Isp*g0));

    rk=rk+Ts*(f1k);
    thetak=thetak+Ts*(f2k);
    phik=phik+Ts*(f3k);
    vrk=vrk+Ts*(f4k);
    vthetak=vthetak+Ts*(f5k);
    vphik=vphik+Ts*(f6k);
    mk=mk+Ts*(f7k);

    xk=[rk thetak phik vrk vthetak vphik mk]';
    
    %%%%%%%%%%%%%%%%%%%% Início da Fase de Predição %%%%%%%%%%%%%%%%%%%%%
    % 1) Predição do Estado a priori utilizando o modelo não linear %%%
%     f1ke=vrke;
%     f2ke=vthetake/(rke*cos(phike));
%     f3ke=vphike/rke;
%     f4ke=((-Tmax*Kk/mke)*sin(betak))+(-um/(rke^2))+(vphike^2/rke)+(vthetake^2/rke)+((rke*w^2)*(cos(phike)))+(2*w*vthetake*cos(phike));
%     f5ke=(((Tmax*Kk)/mke)*cos(betak)*cos(alphak))-(vthetake*vrke/rke)+((vphike*vthetake/rke)*sin(phike)/cos(phike))+(2*w*vphike*sin(phike))+(-2*w*vrke*cos(phike));
%     f6ke=(((Tmax*Kk)/mke)*cos(betak)*sin(alphak))+(-vphike*vrke/rke)+(((-vthetake^2)*sin(phike))/(cos(phike)*rke))+((-rke*w^2)*sin(phike)*cos(phike))+(-2*w*vthetake*sin(phike));
%     f7ke=(-Tmax*Kk/(Isp*g0));
%     
%     rke=rke+Ts*(f1ke);
%     thetake=thetake+Ts*(f2ke);
%     phike=phike+Ts*(f3ke);
%     vrke=vrke+Ts*(f4ke);
%     vthetake=vthetake+Ts*(f5ke);
%     vphike=vphike+Ts*(f6ke);
%     mke=mke+Ts*(f7ke);

    
    % 1) Predição do Estado a priori utilizando o modelo linearizado %%%
    xke=[rke thetake phike vrke vthetake vphike mke]';
    uk=[alphak betak Kk]';
    xke=Ak*xke+Bk*uk;

    % 2) Predição da Covariância %%%%%%
     Pk=Ak*Pk*Ak' + Qk;

    %%%%%%%%%%%%%%%%%%%% Final da Fase de Predição %%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%% Início da Fase de Atualização/Correção %%%%%%%%%%
    
    % 1) Resíduo de Medição %%%%%%%%%%%
    erroy=C*(xk-xke);

    % 2) Resíduo da Covariância %%%%%%%
    Sk=C*Pk*C' + Rk;

    % 3) Ganho Ótimo de Kalman %%%
    Lk=Pk*C'/(Sk);

    % 4) Estado atualizado/corrigido %%%
    xke=xke+Lk*erroy;

    % 5) Covariância Estimada a posteriori %%%%
    Pk=(eye(7)-Lk*C)*Pk*(eye(7)-Lk*C)' + Lk*Rk*Lk';

    %%%%%%%%%%% Final da Fase de Atualização %%%%%%%%%%

    %%%%%%%%%%% Lei de Controle %%%%%%%%%%

      %uk=C2*ref-KLQR*(xke);
      %uk=(-(-KLQR*xke-C2*ref)/(1.157e6)); % 1.157e6@rf=1737,38 e 1.1567e6@rf=1737,52
      %uk=(-(KLQR*xke+C2*ref)/(1.10396e6)); % 1.10396e6@rf=1737,40097 e 1.104e6@rf=1737,428
      %uk=(-(KLQR*(xke+(ref-C*xk)))/(3.1159625e7));
      %uk=(-(KLQR*xke+C2*ref)/(1.00159625e6));
      %uk=(-(KLQR*xke+C2*ref)/(0.90356e6));
      
    %%% Abaixo temos os melhores resultados para a Lei de Controle %%%%
    Kajust1=10.1;%10.03188;%10.01785;%10.0025;%9.995;
    Kajust2=0.986675702480599585; %Kajust2=1;
    Kajust3=2.902e6;
    uk=-(Kajust1*KLQR*xke-Kajust2*C2*(ref-C*rk))/(Kajust3); % Kajust=10.0188; rf=1737.4000235939420235808938741684
%   Kajust1=5.807755217135;%5.775;%5.80755217135;
%   Kajust2=5.809;
%   uk=(-(Kajust1*KLQR*xke+Kajust2*C2*(ref-C*rk))/(2.902065e6)); % Kajust=5.80755217135; rf=1737.4000522598496445425553247333
%     
    alphak=uk(1,1);
    betak=uk(2,1);
    Kk=uk(3,1);
    
    rk=xk(1,1); thetak=xk(2,1); phik=xk(3,1); vrk=xk(4,1); vthetak=xk(5,1); vphik=xk(6,1); mk=xk(7,1);
    rke=xke(1,1); thetake=xke(2,1); phike=xke(3,1); vrke=xke(4,1); vthetake=xke(5,1); vphike=xke(6,1); mke=xke(7,1);

   
    % Armazenando os resultados %%%
    R=[R rk];THETA=[THETA thetak];PHI=[PHI phik];
    VR=[VR vrk];VTHETA=[VTHETA vthetak]; VPHI=[VPHI vphik];M=[M mk];
    Xe=[Xe xke]; W=[W w];

    %if vrk<100 && vrk>-100
        %RKVK=[RKVK rkvk];
    
    %end
end

%%% Plot dos Estados %%%%%%%%

figure
plot(T,R/1000,'m',T_PLOT,Y_PLOT(:,1),'b')
xlabel('t (s)')
ylabel('r (km)')
legend('Distância radial com LQG','Distância radial em malha aberta')
%title('Posição Radial - LQG')
%grid

% figure
% plot(TREF,REF(1,:)/1000,'r')
% xlabel('t (s)')
% ylabel('r em Km')
% title('Posição Radial da Referência - LQG')
 
% figure
% plot(T,((R/1000)-1737.4),'m')
% xlabel('t (s)')
% ylabel('Altitude (km)')
% title('Altitude - LQG')
% %grid

figure
plot(T,THETA*180/pi,'m',T_PLOT,Y_PLOT(:,2)*180/pi,'b')
xlabel('t (s)')
ylabel('theta (°)')
legend('Longitude com LQG','Longitude em malha aberta')
%title('Longitude - LQG')
%grid

figure
plot(T,PHI*180/pi,'m',T_PLOT,Y_PLOT(:,3)*180/pi,'b')
xlabel('t (s)')
ylabel('phi (°)')
legend('Latitude com LQG','Latitude em malha aberta')
%title('Latitude - LQG')
%grid

figure
plot(T,VR,'m',T_PLOT,Y_PLOT(:,4)*1000,'b')
xlabel('t (s)')
ylabel('vr (m/s)')
legend('Velocidade radial com LQG','Velocidade radial em malha aberta')
%title('Velocidade Radial - LQG')
%grid

figure
plot(T,VTHETA,'m',T_PLOT,Y_PLOT(:,5)*1000,'b')
xlabel('t (s)')
ylabel('vtheta (m/s)')
legend('Velocidade longitudinal com LQG','Velocidade longitudinal em malha aberta')
%title('Velocidade Longitudinal - LQG')
%grid

figure
plot(T,VPHI,'m',T_PLOT,Y_PLOT(:,6)*1000,'b')
xlabel('t (s)')
ylabel('vphi (m/s)')
legend('Velocidade latitudinal com LQG','Velocidade latitudinal em malha aberta')
%title('Velocidade Latitudinal - LQG')
%grid

figure
plot(T,M,'m',T_PLOT,Y_PLOT(:,7),'b')
xlabel('t (s)')
ylabel('m (kg)')
legend('Fluxo de massa com LQG','Fluxo de massa em malha aberta')
%title('Massa do Módulo - LQG')
%grid

figure
plot(THETA*180/pi,PHI*180/pi,'m',Y_PLOT(:,2)*180/pi,Y_PLOT(:,3)*180/pi,'b')
xlabel('Longitude em graus')
ylabel('Latitude em graus')
legend('LQG','Malha aberta')
title('Gráfico das Fases - Latitude x Longitude')
% %grid
 
figure
plot3(PHI*180/pi, -THETA*180/pi, (R/1000)-1737.4,'m',Y_PLOT(:,3)*180/pi,-Y_PLOT(:,2)*180/pi,Y_PLOT(:,1)-1737.4,'b')
xlabel('Longitude em graus')
ylabel('Latitude em graus')
zlabel('Altitude em Km')
legend('LQG','Malha aberta')
title('Gráfico das Fases - Altitude x Latitude x Longitude')
% %grid

%vpa([(min(R/1000)); (R(1,end)/1000)])
