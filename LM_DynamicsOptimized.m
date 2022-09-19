%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Projeto Final do curso de Engenharia de Controle e Automação
% Universidade: CEFET - RJ/Uned NI
% Aluna: Laís Lima - Matrícula: 1620368ECAN
% Professor orientador: Mauro Vasconcellos
% Referência principal: Artigo "Three-Dimensional Trajectory Optimization of Soft Lunar Landings from the Parking Orbit with Considerations of the Landing Site" escrito por Bong-Gyun Park and Min-Jea Tahk (2011)
% Script: LM_DynamicsOptimized.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Este código deve ser rodado em primeira instância.

%% Seção 1: Análise em malha aberta
% Nesta seção do código, buscam-se pontos de operação para linearização do sistema 3D de um Módulo Lunar, através de aproximações das curvas obtidas no artigo de referência, cuja abordagem é um controle ótimo em malha aberta.

% % Limpando o espaço de trabalho, a janela de comando e área de figuras do programa
% clc; close all; clear;
% 
% % Parâmetros iniciais e finais para a simulação do sistema dinâmico (resultantes das limitações impostas de rfinal=1737.4km e vrfinal=~+/-1m/s)
% Tff=3400; Tf=3350; Tf2=3340; Tf1=2926; Tf0=2862; rr0=1837.4; %Rf=1737.3979609789;  Vvrf=0,844185604183028 m/s; 
% 
% % Intervalos temporais referentes a cada fase
% span1=linspace(0,7.37,7.37); %old: [0 7.37]
% span2=linspace(7.37, 7.37+Tf0,Tf0); %old: [7.37 7.37+Tf0]
% spanp1=linspace(7.37+Tf0,Tf1,Tf1-Tf0-7.37); %old: [7.37+Tf0 Tf1]
% spanp2=linspace(Tf1,Tf2,Tf2-Tf1); %old: [Tf1 Tf2]
% spanp3=linspace(Tf2,Tf,Tf-Tf2); %old: [Tf2 Tf]
% spanp4=linspace(Tf,Tff,Tff-Tf); %old: [Tf Tff]
% 
% % Representações do sistema dinâmico do módulo lunar nas três fases do problema de aterrissagem
% % De-Orbit Phase (Fase de queima de saída de órbita)
% [t1,y1] = ode45(@(t,x)DeOrbitFunction(t,x,0,(-205+(25/7.37)*t)*(pi/180),1),span1,[rr0;-142.066*(pi/180);-20.24*(pi/180);0;1627.8e-3;-60.139e-3;600]);
%                                                                                  
% % Transfer Phase (Fase de transferência de órbita)
% [t2,y2] = ode45(@(t,x)TransferOrbitFunction(t,x,0,-180*(pi/180),0),span2,y1(end,:)');
% 
% % Powered Descent Phase (Fase de descida motorizada) - Terceira abordagem
% % Esta fase foi subdividida em três partes para melhor ajustar a aproximação das entradas (alpha, betha, k)
% [t3,y3] = ode45(@(t,x)PDFp1(t,x,Tf0,Tf1,Tf),spanp1,y2(end,:)'); % Powered Descent Phase - part 1                                                             
% [t4,y4] = ode45(@(t,x)PDFp2(t,x,Tf0,Tf1,Tf2,Tf),spanp2,y3(end,:)'); % Powered Descent Phase - part 2                                                              
% [t5,y5] = ode45(@(t,x)PDFp3(t,x,Tf0,Tf2,Tf),spanp3,y4(end,:)'); % Powered Descent Phase - part 3                                                             
% [t6,y6] = ode45(@(t,x)PDFp4(t,x),spanp4,y5(end,:)');
% % Reunindo as respostas do sistema
% t=[t1;t2;t3;t4;t5;t6];
% y=[y1;y2;y3;y4;y5;y6];

%% Seção 2: Linearização
% Após a simulação em malha aberta e a reunião de pontos de operação, faz-se necessário linearizar o sistema 3D por partes, isto é, em cada ponto de operação, de forma que o projeto de controladores ótimos LQR e LQG seja possível.

% Criando as variáveis simbólicas a serem utilizadas
syms r theta phi vr vtheta vphi m alpha betha k x1 x2 x3 x4 x5 x6 x7 u1 u2 u3 Tmax mi omega Isp g0

% Vetor de estados
x = [r theta phi vr vtheta vphi m]'; %x=[x1 x2 x3 x4 x5 x6 x7]';
r = x1; % Distância radial
theta = x2; % Longitude
phi = x3; % Latitude
vr = x4; % Velocidade radial/vertical
vtheta = x5; % Velocidade Longitudinal
vphi = x6; % Velocidade Latitudinal
m = x7; % Massa

% Vetor de entradas
u = [alpha betha k]'; %=[u1 u2 u3]';
alpha = u1; % Ângulo de direção do propulsor no plano local
betha = u2; % Ângulo de direção do propulsor acima do plano local
k = u3; % Acelerador/propulsor

% Equações dinâmicas
f1 = vr; %r_ponto
f2 = vtheta/(r*cos(phi)); %theta_ponto
f3 = vphi/r; %phi_ponto
f4 = ((-Tmax*k)/m)*sin(betha) - (mi)/(r^2) + ((vphi^2))/r + ((vtheta)^2)/r + r*(omega^2)*cos(phi) + 2*omega*vtheta*cos(phi); %vr_ponto
f5 = ((Tmax*k)/m)*cos(betha)*cos(alpha) - (vr*vtheta)/r + (vtheta*vphi*sin(phi))/(r*cos(phi)) + 2*omega*vphi*sin(phi) - 2*omega*vr*cos(phi); %vtheta_ponto
f6 = ((Tmax*k)/m)*cos(betha)*sin(alpha) - (vr*vphi)/r - ((vtheta^2)*sin(phi))/(r*cos(phi)) - r*(omega^2)*sin(phi)*cos(phi) - 2*omega*vtheta*sin(phi); %vphi_ponto
f7 = -(Tmax*k)/(Isp*g0); %m_ponto

% Cálculo dos jacobianos para x e u
Jx=[diff(f1,x1) diff(f1,x2) diff(f1,x3) diff(f1,x4) diff(f1,x5) diff(f1,x6) diff(f1,x7);
    diff(f2,x1) diff(f2,x2) diff(f2,x3) diff(f2,x4) diff(f2,x5) diff(f2,x6) diff(f2,x7);
    diff(f3,x1) diff(f3,x2) diff(f3,x3) diff(f3,x4) diff(f3,x5) diff(f3,x6) diff(f3,x7);
    diff(f4,x1) diff(f4,x2) diff(f4,x3) diff(f4,x4) diff(f4,x5) diff(f4,x6) diff(f4,x7);
    diff(f5,x1) diff(f5,x2) diff(f5,x3) diff(f5,x4) diff(f5,x5) diff(f5,x6) diff(f5,x7);
    diff(f6,x1) diff(f6,x2) diff(f6,x3) diff(f6,x4) diff(f6,x5) diff(f6,x6) diff(f6,x7);
    diff(f7,x1) diff(f7,x2) diff(f7,x3) diff(f7,x4) diff(f7,x5) diff(f7,x6) diff(f7,x7)];

Ju=[diff(f1,u1) diff(f1,u2) diff(f1,u3);
    diff(f2,u1) diff(f2,u2) diff(f2,u3);
    diff(f3,u1) diff(f3,u2) diff(f3,u3);
    diff(f4,u1) diff(f4,u2) diff(f4,u3);
    diff(f5,u1) diff(f5,u2) diff(f5,u3);
    diff(f6,u1) diff(f6,u2) diff(f6,u3);
    diff(f7,u1) diff(f7,u2) diff(f7,u3)];

% Inicializando vetores que armazenarão os parâmetros para a montagem dos gráficos
F1=[]; % Representa a primeira equação dinâmica
F2=[]; % Representa a segunda equação dinâmica
F3=[]; % Representa a terceira equação dinâmica
F4=[]; % Representa a quarta equação dinâmica
F5=[]; % Representa a quinta equação dinâmica
F6=[]; % Representa a sexta equação dinâmica
F7=[]; % Representa a sétima equação dinâmica
A_lm=[]; % Representa a matriz de estados A do sistema linearizado
B_lm=[]; % Representa a matriz de entradas B do sistema linearizado
T=[]; % Vetor de instantes
R_LM=[]; % Vetor do estado r
THETA_LM=[]; % Vetor do estado theta
PHI_LM=[]; % Vetor do estado phi
VR_LM=[]; % Vetor do estado vr
VTHETA_LM=[]; % Vetor do estado vtheta
VPHI_LM=[]; % Vetor do estado vphi
M_LM=[]; % Vetor do estado m
ALPHA_LM=[]; % Vetor da entrada alpha
BETHA_LM=[]; % Vetor da entrada betha
K_LM=[]; % Vetor da entrada k

% Matrizes de parâmetros controláveis
T_valido=[]; 
R_valido=[]; 
THETA_valido=[];
PHI_valido=[]; 
VR_valido=[]; 
VTHETA_valido=[]; 
VPHI_valido=[]; 
M_valido=[]; 
ALPHA_valido=[];
BETHA_valido=[];
K_valido=[];
A_valida=[];
B_valida=[];

% Matrizes de parâmetros observáveis
T_observavel=[];
R_observavel=[];
THETA_observavel=[];
PHI_observavel=[]; 
VR_observavel=[]; 
VTHETA_observavel=[]; 
VPHI_observavel=[]; 
M_observavel=[]; 
ALPHA_observavel=[];
BETHA_observavel=[];
K_observavel=[];
A_observavel=[];
B_observavel=[];

% Preenchendo os vetores de parâmetros calculados
for i=1:1:length(t)
     % Selecionando os pontos de operação atrelados aos seus respectivos instantes de tempo
     r_lm=y(i,1); % [r]=km
     theta_lm=y(i,2); % [theta]=rad
     phi_lm=y(i,3); % [phi]=rad
     vr_lm=y(i,4); % [vr]=km/s
     vtheta_lm=y(i,5); % [vtheta]=km/s
     vphi_lm=y(i,6); % [vphi]=km/s
     m_lm=y(i,7); % [m]=kg
     
     if (t(i)>=0 && t(i)<=7.37)
         alpha_lm=0; % [alpha]=rad
         betha_lm=(-205+(25/7.37)*t(i))*(pi/180); % [betha]=rad
         k_lm=1; % adimensional
     end
     if (t(i)>7.37 && t(i)<=7.37+Tf0)
         alpha_lm=0; % [alpha]=rad
         betha_lm=-180*(pi/180); % [betha]=rad
         k_lm=0; %modificar para 10^(-5), verificar 10^(-6), 10^(-4) % adimensional
         %k_lm=1e-5; % adimensional
     end
     if (t(i)>7.37+Tf0) && (t(i)<=Tf1) 
         betha_lm=(-180+30*(t(i)-Tf0)/(Tf-Tf0))*pi/180; % [betha]=rad
         alpha_lm=(6.5*(t(i)-Tf0)/(Tf1-Tf0))*pi/180; % [alpha]=rad
         k_lm=1; % adimensional
     end
     if (t(i)>Tf1) && (t(i)<=Tf2)
         betha_lm=(-180+30*(t(i)-Tf0)/(Tf-Tf0))*pi/180; % [betha]=rad
         alpha_lm=(6.5-0.5*(t(i)-Tf1)/(Tf2-Tf1))*pi/180; % [alpha]=rad
         k_lm=1; % adimensional
     end
     if (t(i)>Tf2) && (t(i)<=Tf)
         betha_lm=(-180+30*(t(i)-Tf0)/(Tf-Tf0))*pi/180; % [betha]=rad
         alpha_lm=(6+3*(t(i)-Tf2)/(Tf-Tf2))*pi/180; % [alpha]=rad
         k_lm=1; % adimensional
     end
     if (t(i)>Tf) && (t(i)<=Tff)
         betha_lm=0; % [betha]=rad
         alpha_lm=0; % [alpha]=rad
         k_lm=0; % adimensional
         %k_lm=1e-5; % adimensional
     end

     % Montando os vetores de estados
     R_LM=[R_LM; r_lm];
     THETA_LM=[THETA_LM; theta_lm];
     PHI_LM=[PHI_LM; phi_lm];
     VR_LM=[VR_LM; vr_lm];
     VTHETA_LM=[VTHETA_LM; vtheta_lm];
     VPHI_LM=[VPHI_LM; vphi_lm];
     M_LM=[M_LM; m_lm];
     % Montando o vetor de instantes
     T=[T t(i)];
     % Montando os vetores de entradas
     ALPHA_LM=[ALPHA_LM; alpha_lm];
     BETHA_LM=[BETHA_LM; betha_lm];
     K_LM=[K_LM; k_lm];
end

% Aplicando a linearização em torno dos pontos de operação
for i=1:1:length(y)
     % Montando as matrizes A e B através dos jacobianos
     A = subs(Jx,[mi,Tmax,omega,Isp,g0,x1,x2,x3,x4,x5,x6,x7,u1,u2,u3],[4902.78,1700*10^(-3),2.6632*10^(-6),316,9.81*10^(-3),R_LM(i),THETA_LM(i),PHI_LM(i),VR_LM(i),VTHETA_LM(i),VPHI_LM(i),M_LM(i),ALPHA_LM(i),BETHA_LM(i),K_LM(i)]);
     B = subs(Ju,[mi,Tmax,omega,Isp,g0,x1,x2,x3,x4,x5,x6,x7,u1,u2,u3],[4902.78,1700*10^(-3),2.6632*10^(-6),316,9.81*10^(-3),R_LM(i),THETA_LM(i),PHI_LM(i),VR_LM(i),VTHETA_LM(i),VPHI_LM(i),M_LM(i),ALPHA_LM(i),BETHA_LM(i),K_LM(i)]);
     A = double(A);
     B = double(B);
     % Verificando a controlabilidade do sistema
     controlabilidade=ctrb(A,B); 
     postoCon=rank(controlabilidade);
     tamanho=size(controlabilidade);
   if(postoCon==tamanho(1)) %comentar esse if para pegar todos os pontos (controláveis e incontroláveis)
     % Reunindo valores válidos que se enquadram no critério de controlabilidade
        T_valido=[T_valido; t(i)];
        R_valido=[R_valido; R_LM(i)];
        THETA_valido=[THETA_valido; THETA_LM(i)];
        PHI_valido=[PHI_valido; PHI_LM(i)]; 
        VR_valido=[VR_valido; VR_LM(i)]; 
        VTHETA_valido=[VTHETA_valido; VTHETA_LM(i)]; 
        VPHI_valido=[VPHI_valido; VPHI_LM(i)];
        M_valido=[M_valido; M_LM(i)]; 
        ALPHA_valido=[ALPHA_valido; ALPHA_LM(i)];
        BETHA_valido=[BETHA_valido; BETHA_LM(i)];
        K_valido=[K_valido; K_LM(i)];
        A_valida=[A_valida; A];
        B_valida=[B_valida; B];
   end
     % Verificando a observabilidade dos sistema
     C=eye(size(A,1),size(A,2));
     observabilidade=obsv(A,C);
     postoObs=rank(observabilidade);
     if(postoObs==tamanho(1))
        T_observavel=[T_observavel; t(i)];
        R_observavel=[R_observavel; R_LM(i)];
        THETA_observavel=[THETA_observavel; THETA_LM(i)];
        PHI_observavel=[PHI_observavel; PHI_LM(i)]; 
        VR_observavel=[VR_observavel; VR_LM(i)]; 
        VTHETA_observavel=[VTHETA_observavel; VTHETA_LM(i)]; 
        VPHI_observavel=[VPHI_observavel; VPHI_LM(i)];
        M_observavel=[M_observavel; M_LM(i)]; 
        ALPHA_observavel=[ALPHA_observavel; ALPHA_LM(i)];
        BETHA_observavel=[BETHA_observavel; BETHA_LM(i)];
        K_observavel=[K_observavel; K_LM(i)];
        A_observavel=[A_observavel; A];
        B_observavel=[B_observavel; B];    
     end
    
     % Montando o conjunto de matrizes A e B do sistema linearizado
     A_lm = [A_lm; A;];
     B_lm = [B_lm; B;];
%      % Calculando as variações do sistema nos pontos de operação
%      f1_calc = subs(f1,vr,VR_LM(i));
%      f2_calc = subs(f2,[vtheta,r,phi],[VTHETA_LM(i),R_LM(i),PHI_LM(i)]);
%      f3_calc = subs(f3,[vphi,r],[VPHI_LM(i),R_LM(i)]);
%      f4_calc = subs(f4,[Tmax,k,m,betha,mi,r,vphi,vtheta,omega,phi],[1700e-3,K_LM(i),M_LM(i),BETHA_LM(i),4902.78,R_LM(i),VPHI_LM(i),VTHETA_LM(i),2.6632*10^(-6),PHI_LM(i)]);
%      f5_calc = subs(f5,[Tmax,k,m,betha,alpha,vr,vtheta,r,vphi,phi,omega],[1700e-3,K_LM(i),M_LM(i),BETHA_LM(i),ALPHA_LM(i),VR_LM(i),VTHETA_LM(i),R_LM(i),VPHI_LM(i),PHI_LM(i),2.6632*10^(-6)]);
%      f6_calc = subs(f6,[Tmax,k,m,betha,alpha,vr,vtheta,r,vphi,phi,omega],[1700e-3,K_LM(i),M_LM(i),BETHA_LM(i),ALPHA_LM(i),VR_LM(i),VTHETA_LM(i),R_LM(i),VPHI_LM(i),PHI_LM(i),2.6632*10^(-6)]);
%      f7_calc = subs(f7,[Tmax,k,Isp,g0],[1700e-3,K_LM(i),316,9.81e-3]);
%      % Montando os vetores de variações do sistema nos pontos de operação
%      F1=[F1 f1_calc];
%      F2=[F2 f2_calc];
%      F3=[F3 f3_calc];
%      F4=[F4 f4_calc];
%      F5=[F5 f5_calc];
%      F6=[F6 f6_calc];
%      F7=[F7 f7_calc];
end

% Matrizes de parâmetros controláveis para PLOT
T_valido_PLOT=[]; 
R_valido_PLOT=[]; 
THETA_valido_PLOT=[];
PHI_valido_PLOT=[]; 
VR_valido_PLOT=[]; 
VTHETA_valido_PLOT=[]; 
VPHI_valido_PLOT=[]; 
M_valido_PLOT=[]; 

% Matrizes de parâmetros observáveis para PLOT
T_observavel_PLOT=[];
R_observavel_PLOT=[];
THETA_observavel_PLOT=[];
PHI_observavel_PLOT=[]; 
VR_observavel_PLOT=[]; 
VTHETA_observavel_PLOT=[]; 
VPHI_observavel_PLOT=[]; 
M_observavel_PLOT=[]; 

for i=1:1:length(T_valido)
    if(T_valido(i)<=tempodeparada)
        T_valido_PLOT=[T_valido_PLOT; T_valido(i)];
        R_valido_PLOT=[R_valido_PLOT; R_valido(i)]; 
        THETA_valido_PLOT=[THETA_valido_PLOT; THETA_valido(i)];
        PHI_valido_PLOT=[PHI_valido_PLOT; PHI_valido(i)]; 
        VR_valido_PLOT=[VR_valido_PLOT; VR_valido(i)]; 
        VTHETA_valido_PLOT=[VTHETA_valido_PLOT; VTHETA_valido(i)]; 
        VPHI_valido_PLOT=[VPHI_valido_PLOT; VPHI_valido(i)];
        M_valido_PLOT=[M_valido_PLOT; M_valido(i)];
    end
end

for i=1:1:length(T_observavel)
    if(T_observavel(i)<=tempodeparada)
        T_observavel_PLOT=[T_observavel_PLOT; T_observavel(i)];
        R_observavel_PLOT=[R_observavel_PLOT; R_observavel(i)]; 
        THETA_observavel_PLOT=[THETA_observavel_PLOT; THETA_observavel(i)];
        PHI_observavel_PLOT=[PHI_observavel_PLOT; PHI_observavel(i)]; 
        VR_observavel_PLOT=[VR_observavel_PLOT; VR_observavel(i)]; 
        VTHETA_observavel_PLOT=[VTHETA_observavel_PLOT; VTHETA_observavel(i)]; 
        VPHI_observavel_PLOT=[VPHI_observavel_PLOT; VPHI_observavel(i)];
        M_observavel_PLOT=[M_observavel_PLOT; M_observavel(i)];
    end
end


%% Seção 3: Gráficos dos estados do sistema na solução em malha aberta

% Distância radial
figure
plot(T_PLOT,Y_PLOT(:,1),'b-')
title('Estado 1 - Distância radial')
ylabel('r (km)')
xlabel('t(s)')
figure
subplot(2,1,1)
plot(T_valido_PLOT,R_valido_PLOT,'o')
title('Estado 1 - Controlabilidade')
ylabel('r (km)')
xlabel('t(s)')
subplot(2,1,2)
plot(T_observavel_PLOT,R_observavel_PLOT,'+')
title('Estado 1 - Observabilidade')
ylabel('r (km)')
xlabel('t(s)')

% Longitude
figure
plot(T_PLOT,Y_PLOT(:,2)*180/pi,'r-')
title('Estado 2 - Longitude')
ylabel('theta (°)')
xlabel('t(s)')
figure
subplot(2,1,1)
plot(T_valido_PLOT,THETA_valido_PLOT*180/pi,'o')
title('Estado 2 - Controlabilidade')
ylabel('theta (°)')
xlabel('t(s)')
subplot(2,1,2)
plot(T_observavel_PLOT,THETA_observavel_PLOT*180/pi,'+')
title('Estado 2 - Observabilidade')
ylabel('theta (°)')
xlabel('t(s)')

% Latitude
figure
plot(T_PLOT,Y_PLOT(:,3)*180/pi,'g-')
title('Estado 3 - Latitude')
ylabel('phi (°)')
xlabel('t(s)')
figure
subplot(2,1,1)
plot(T_valido_PLOT,PHI_valido_PLOT*180/pi,'o')
title('Estado 3 - Controlabilidade')
ylabel('phi (°)')
xlabel('t(s)')
subplot(2,1,2)
plot(T_observavel_PLOT,PHI_observavel_PLOT*180/pi,'+')
title('Estado 3 - Observabilidade')
ylabel('phi (°)')
xlabel('t(s)')

% Velocidade radial
figure
plot(T_PLOT,Y_PLOT(:,4)*1000,'m-')
title('Estado 4 - Veloc. Radial')
ylabel('vr (m/s)')
xlabel('t(s)')
figure
subplot(2,1,1)
plot(T_valido_PLOT,VR_valido_PLOT*1000,'o')
title('Estado 4 - Controlabilidade')
ylabel('vr (m/s)')
xlabel('t(s)')
subplot(2,1,2)
plot(T_observavel_PLOT,VR_observavel_PLOT*1000,'+')
title('Estado 4 - Observabilidade')
ylabel('vr (m/s)')
xlabel('t(s)')

% Velocidade longitudinal
figure
plot(T_PLOT,Y_PLOT(:,5)*1000,'y-')
title('Estado 5 - Veloc. Longitudinal')
ylabel('vtheta (m/s)')
xlabel('t(s)')
figure
subplot(2,1,1)
plot(T_valido_PLOT,VTHETA_valido_PLOT*1000,'o')
title('Estado 5 - Controlabilidade')
ylabel('vtheta (m/s)')
xlabel('t(s)')
subplot(2,1,2)
plot(T_observavel_PLOT,VTHETA_observavel_PLOT*1000,'+')
title('Estado 5 - Observabilidade')
ylabel('vtheta (m/s)')
xlabel('t(s)')

% Velocidade latitudinal
figure
plot(T_PLOT,Y_PLOT(:,6)*1000,'c-')
title('Estado 6 - Veloc. Latitudinal')
ylabel('vphi (m/s)')
xlabel('t(s)')
figure
subplot(2,1,1)
plot(T_valido_PLOT,VPHI_valido_PLOT*1000,'o')
title('Estado 6 - Controlabilidade')
ylabel('vphi (m/s)')
xlabel('t(s)')
subplot(2,1,2)
plot(T_observavel_PLOT,VPHI_observavel_PLOT*1000,'+')
title('Estado 6 - Observabilidade')
ylabel('vphi (m/s)')
xlabel('t(s)')

% Massa
figure
plot(T_PLOT,Y_PLOT(:,7),'k-')
title('Estado 7 - Massa')
ylabel('m (kg)')
xlabel('t(s)')
figure
subplot(2,1,1)
plot(T_valido_PLOT,M_valido_PLOT,'o')
title('Estado 7 - Controlabilidade')
ylabel('m (kg)')
xlabel('t(s)')
subplot(2,1,2)
plot(T_observavel_PLOT,M_observavel_PLOT,'+')
title('Estado 7 - Observabilidade')
ylabel('m (kg)')
xlabel('t(s)')

% Gráfico de fase (Latitude x Longitude)
figure
plot(Y_PLOT(:,2)*180/pi,Y_PLOT(:,3)*180/pi,'m-')
title('Gráfico de Fase')
xlabel('Longitude (°)')
ylabel('Latitude (°)')

% Altitude x Longitude x Latitude
figure
plot3(Y_PLOT(:,3)*180/pi,Y_PLOT(:,2)*180/pi,Y_PLOT(:,1)-1737.4)
xlabel('Latitude (°)')
ylabel('Longitude (°)')
zlabel('Altitude (km)')
title('Gráfico Altitude x Longitude x Latitude')


%% Seção 3: Gráfico das variações dos estados da solução em malha aberta
% figure
% plot(T,F1,'b-*')
% ylabel('rPonto(km)')
% xlabel('tempo(s)')
% title('Estado 1 - Variação da Posição Radial')
% 
% figure
% plot(T,F2*180/pi,'r-*')
% ylabel('thetaPonto(°)')
% xlabel('tempo(s)')
% title('Estado 2 - Variação da Longitude')
% 
% figure
% plot(T,F3*180/pi,'g-*')
% ylabel('phiPonto(°)')
% xlabel('tempo(s)')
% title('Estado 3 - Variação da Latitude')
% 
% figure
% plot(T,F4*1000,'y-*')
% ylabel('vrPonto(m/s)')
% xlabel('tempo(s)')
% title('Estado 4 - Aceleração Radial')
% 
% figure
% plot(T,F5*1000,'m-*')
% ylabel('vthetaPonto(m/s)')
% xlabel('tempo(s)')
% title('Estado 5 - Aceleração Longitudinal')
% 
% figure
% plot(T,F6*1000,'c-*')
% ylabel('vphiPonto(m/s)')
% xlabel('tempo(s)')
% title('Estado 6 - Aceleração Latitudinal')
% 
% figure
% plot(T,F7,'k-*')
% ylabel('mPonto(kg)')
% xlabel('tempo(s)')
% title('Estado 7 - Fluxo de Massa')

% %% Seção 4: Funções representativas da dinâmica do sistema 3D em cada fase do problema de aterrissagem
% % De-Orbit Phase
% function dydt = DeOrbitFunction(t,x,alpha,betha,k)
% g0=9.81e-3; % Gravidade da Terra em Km/s^2
% Isp=316; % Impulso Específico em segundos
% Tmax=1700e-3; % Força de Propulsão Máxima em Kg*Km/s^2
% mi=4902.78; % Parâmetro Gravitacional Padrão em km^3/s^2
% omega=2.6632e-6; % Velocidade Angular da Lua em rad/s.
% dydt=[x(4); 
%       x(5)/(x(1)*cos(x(3))); 
%       x(6)/x(1); 
%       ((-Tmax*k)/x(7))*sin(betha) - (mi)/(x(1)^2) + (x(6)^2)/x(1) + (x(5)^2)/x(1) + x(1)*(omega^2)*cos(x(3)) + 2*omega*x(5)*cos(x(3));
%       ((Tmax*k)/x(7))*cos(betha)*cos(alpha) - (x(4)*x(5))/x(1) + (x(5)*x(6)*sin(x(3)))/(x(1)*cos(x(3))) + 2*omega*x(6)*sin(x(3)) - 2*omega*x(4)*cos(x(3));
%       ((Tmax*k)/x(7))*cos(betha)*sin(alpha) - (x(4)*x(6))/x(1) - ((x(5)^2)*sin(x(3)))/(x(1)*cos(x(3))) - x(1)*(omega^2)*sin(x(3))*cos(x(3)) - 2*omega*x(5)*sin(x(3));
%       -(Tmax*k)/(Isp*g0)];
% end
% 
% % Transfer Orbit Phase
% function dydt = TransferOrbitFunction(t,x,alpha,betha,k)
% g0=9.81e-3; % Gravidade da Terra em Km/s^2
% Isp=316; % Impulso Específico em segundos
% Tmax=1700e-3; % Força de Propulsão Máxima em Kg*Km/s^2
% mi=4902.78; % Parâmetro Gravitacional Padrão em km^3/s^2
% omega=2.6632e-6; % Velocidade Angular da Lua em rad/s.
% dydt=[x(4); 
%       x(5)/(x(1)*cos(x(3))); 
%       x(6)/x(1); 
%       ((-Tmax*k)/x(7))*sin(betha) - (mi)/(x(1)^2) + ((x(6)^2))/x(1) + ((x(5))^2)/x(1) + x(1)*(omega^2)*cos(x(3)) + 2*omega*x(5)*cos(x(3));
%       ((Tmax*k)/x(7))*cos(betha)*cos(alpha) - (x(4)*x(5))/x(1) + (x(5)*x(6)*sin(x(3)))/(x(1)*cos(x(3))) + 2*omega*x(6)*sin(x(3)) - 2*omega*x(4)*cos(x(3));
%       ((Tmax*k)/x(7))*cos(betha)*sin(alpha) - (x(4)*x(6))/x(1) - ((x(5)^2)*sin(x(3)))/(x(1)*cos(x(3))) - x(1)*(omega^2)*sin(x(3))*cos(x(3)) - 2*omega*x(5)*sin(x(3));
%       -(Tmax*k)/(Isp*g0)];
% end
% 
% %%% Powered Descent Phase - Terceira abordagem
% % Powered Descent Phase - part 1
% function dydt = PDFp1(t,x,Tf0,Tf1,Tf)
% betha=(-180+30*(t-Tf0)/(Tf-Tf0))*pi/180;
% alpha=(6.5*(t-Tf0)/(Tf1-Tf0))*pi/180;
% k=1;
% %k=0.0025*t-7.0093; %Aproximei pegando 2 pontos do gráfico (ver bloco de notas das condições iniciais)
% g0=9.81e-3; % gravidade da Terra em Km/s^2
% Isp=316; % Impulso Específico em segundos
% Tmax=1700e-3; % Força de Propulsão Máxima em Kg*Km/s^2
% mi=4902.78; % Parâmetro Gravitacional Padrão em km^3/s^2
% omega=2.6632e-6; % Velocidade Angular da Lua em rad/s.
% dydt=[x(4); 
%       x(5)/(x(1)*cos(x(3))); 
%       x(6)/x(1); 
%       ((-Tmax*k)/x(7))*sin(betha) - (mi)/(x(1)^2) + ((x(6)^2))/x(1) + ((x(5))^2)/x(1) + x(1)*(omega^2)*cos(x(3)) + 2*omega*x(5)*cos(x(3));
%       ((Tmax*k)/x(7))*cos(betha)*cos(alpha) - (x(4)*x(5))/x(1) + (x(5)*x(6)*sin(x(3)))/(x(1)*cos(x(3))) + 2*omega*x(6)*sin(x(3)) - 2*omega*x(4)*cos(x(3));
%       ((Tmax*k)/x(7))*cos(betha)*sin(alpha) - (x(4)*x(6))/x(1) - ((x(5)^2)*sin(x(3)))/(x(1)*cos(x(3))) - x(1)*(omega^2)*sin(x(3))*cos(x(3)) - 2*omega*x(5)*sin(x(3));
%       -(Tmax*k)/(Isp*g0)];
% end
% 
% % Powered Descent Phase - part 2
% function dydt = PDFp2(t,x,Tf0,Tf1,Tf2,Tf)
% betha=(-180+30*(t-Tf0)/(Tf-Tf0))*pi/180;
% alpha=(6.5-0.5*(t-Tf1)/(Tf2-Tf1))*pi/180;
% k=1;
% g0=9.81e-3; % gravidade da Terra em Km/s^2
% Isp=316; % Impulso Específico em segundos
% Tmax=1700e-3; % Força de Propulsão Máxima em Kg*Km/s^2
% mi=4902.78; % Parâmetro Gravitacional Padrão em km^3/s^2
% omega=2.6632e-6; % Velocidade Angular da Lua em rad/s.
% dydt=[x(4); 
%       x(5)/(x(1)*cos(x(3))); 
%       x(6)/x(1); 
%       ((-Tmax*k)/x(7))*sin(betha) - (mi)/(x(1)^2) + ((x(6)^2))/x(1) + ((x(5))^2)/x(1) + x(1)*(omega^2)*cos(x(3)) + 2*omega*x(5)*cos(x(3));
%       ((Tmax*k)/x(7))*cos(betha)*cos(alpha) - (x(4)*x(5))/x(1) + (x(5)*x(6)*sin(x(3)))/(x(1)*cos(x(3))) + 2*omega*x(6)*sin(x(3)) - 2*omega*x(4)*cos(x(3));
%       ((Tmax*k)/x(7))*cos(betha)*sin(alpha) - (x(4)*x(6))/x(1) - ((x(5)^2)*sin(x(3)))/(x(1)*cos(x(3))) - x(1)*(omega^2)*sin(x(3))*cos(x(3)) - 2*omega*x(5)*sin(x(3));
%       -(Tmax*k)/(Isp*g0)];
% end
% 
% % Powered Descent Phase - part 3
% function dydt = PDFp3(t,x,Tf0,Tf2,Tf)
% betha=(-180+30*(t-Tf0)/(Tf-Tf0))*pi/180;
% alpha=(6+3*(t-Tf2)/(Tf-Tf2))*pi/180;
% k=1;
% g0=9.81e-3; % km/s^2
% Isp=316; % segundos
% Tmax=1700e-3; % kg*km/s^2
% mi=4902.78; % km^3/s^2
% omega=2.6632e-6; % rad/s.
% dydt=[x(4); 
%       x(5)/(x(1)*cos(x(3))); 
%       x(6)/x(1); 
%       ((-Tmax*k)/x(7))*sin(betha) - (mi)/(x(1)^2) + ((x(6)^2))/x(1) + ((x(5))^2)/x(1) + x(1)*(omega^2)*cos(x(3)) + 2*omega*x(5)*cos(x(3));
%       ((Tmax*k)/x(7))*cos(betha)*cos(alpha) - (x(4)*x(5))/x(1) + (x(5)*x(6)*sin(x(3)))/(x(1)*cos(x(3))) + 2*omega*x(6)*sin(x(3)) - 2*omega*x(4)*cos(x(3));
%       ((Tmax*k)/x(7))*cos(betha)*sin(alpha) - (x(4)*x(6))/x(1) - ((x(5)^2)*sin(x(3)))/(x(1)*cos(x(3))) - x(1)*(omega^2)*sin(x(3))*cos(x(3)) - 2*omega*x(5)*sin(x(3));
%       -(Tmax*k)/(Isp*g0)];
% end
% 
% % Powered Descent Phase - part 4
% function dydt = PDFp4(t,x)
% betha=0;
% alpha=0;
% k=0;
% g0=9.81e-3; % km/s^2
% Isp=316; % segundos
% Tmax=1700e-3; % kg*km/s^2
% mi=4902.78; % km^3/s^2
% omega=2.6632e-6; % rad/s.
% dydt=[x(4); 
%       x(5)/(x(1)*cos(x(3))); 
%       x(6)/x(1); 
%       ((-Tmax*k)/x(7))*sin(betha) - (mi)/(x(1)^2) + ((x(6)^2))/x(1) + ((x(5))^2)/x(1) + x(1)*(omega^2)*cos(x(3)) + 2*omega*x(5)*cos(x(3));
%       ((Tmax*k)/x(7))*cos(betha)*cos(alpha) - (x(4)*x(5))/x(1) + (x(5)*x(6)*sin(x(3)))/(x(1)*cos(x(3))) + 2*omega*x(6)*sin(x(3)) - 2*omega*x(4)*cos(x(3));
%       ((Tmax*k)/x(7))*cos(betha)*sin(alpha) - (x(4)*x(6))/x(1) - ((x(5)^2)*sin(x(3)))/(x(1)*cos(x(3))) - x(1)*(omega^2)*sin(x(3))*cos(x(3)) - 2*omega*x(5)*sin(x(3));
%       -(Tmax*k)/(Isp*g0)];
% end
