%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Projeto Final do curso de Engenharia de Controle e Automação
% Universidade: CEFET - RJ/Uned NI
% Aluna: Laís Lima - Matrícula: 1620368ECAN
% Professor orientador: Mauro Vasconcellos
% Referência principal: Artigo "Three-Dimensional Trajectory Optimization of Soft Lunar Landings from the Parking Orbit with Considerations of the Landing Site" escrito por Bong-Gyun Park and Min-Jea Tahk (2011)
% Script: LM_LQR1.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Projetando o controlador LQR por tentativa e erro para o sistema 3D de um módulo lunar 

%% Seção 1 - Definindo as matrizes de ponderação

% Limpando a área de gráficos
%close all;

% Definindo pesos componentes da matriz Q e R por tentativa e erro

% Pesos antigos (prof Mauro)
%Q1=1; Q2=1; Q3=1; Q4=1; Q5=1; Q6=1; Q7=1; %não está bom
%Q1=.03; Q2=1; Q3=1; Q4=4.500; Q5=10; Q6=10; Q7=20; %não está bom
%Q1=1; Q2=1; Q3=1; Q4=10; Q5=10; Q6=10; Q7=.5; % não está bom
%Q1=.03; Q2=1; Q3=1; Q4=4.500; Q5=10; Q6=10; Q7=.5; % r, vr muito distantes, m pouco distante na fase de transferência
%Q1=.03; Q2=1; Q3=1; Q4=4500; Q5=10; Q6=10; Q7=20; %não está bom

% Minhas novas tentativas
%Q1=.03; Q2=1; Q3=10; Q4=4.5; Q5=10; Q6=10; Q7=.5; % começou a melhorar todos os estados
%Q1=.03; Q2=1; Q3=13.5; Q4=4.5; Q5=10; Q6=10; Q7=.5; % muito bom (mas vr
%ainda tá ruim no final)
%Q1=.09; Q2=1; Q3=13.5; Q4=5; Q5=10; Q6=10; Q7=.5; % muito muito bom (mas vr
%ainda tá ruim no final)
Q1=.09; Q2=1; Q3=13.5; Q4=0.2; Q5=10; Q6=10; Q7=.5; 

Q=[Q1 0 0 0 0 0 0;
   0 Q2 0 0 0 0 0;
   0 0 Q3 0 0 0 0;
   0 0 0 Q4 0 0 0;
   0 0 0 0 Q5 0 0;
   0 0 0 0 0 Q6 0;
   0 0 0 0 0 0 Q7];
R= diag([1,1,500]);

% Criando variáveis auxiliares
ControleU=[];
y_valido=[R_valido THETA_valido PHI_valido VR_valido VTHETA_valido VPHI_valido M_valido];
y_transposto=y_valido';
estadosX=[];
Matrizes_A=[];
Matrizes_B=[];

ref=[1737.4; 0.3546;0.3575;-1e-3;0;0;250]; %estado de referência/estado final

%Reescrevendo  os vetores de estados em uma matriz (7 x m), em que 
%m=length(y_válido). Apenas os estados válidos do sistema foram
%considerados, de acordo com o critério de controlabilidade

for i=1:1:length(y_valido)
    estadosX=[estadosX; y_transposto(:,i)];
end


%% Seção 2: Calculando as leis de controle 
inc=2;

for i=1:7:(length(A_valida)-(7*inc+6))
        A_atual=A_valida(i:i+6,:);
        B_atual=B_valida(i:i+6,:);
        % Calculando o ganho K
        GanhoK_atual=lqr(A_atual,B_atual,Q,R);
        % Montando a lei de controle LQR
        ref0=estadosX(i+7*inc:i+7*inc+6,:);
        ControleU_atual=-GanhoK_atual*(estadosX(i:i+6,:)-ref0);      
        ControleU=[ControleU ControleU_atual];
end
 
cntdr=0;
tempos_lqr1=[];

%% Seção 3: Aplicando as leis de controle no sistema 3D
%for i=1:inc:(length(ControleU)-(inc+1))
for i=1:inc:(length(ControleU)-(inc+1))
    cntdr=cntdr+1;
    if(T_valido(i)>=0 && T_valido(i)<7.37)
        t1_lqr_span = linspace(T_valido(i),T_valido(i+inc),2);
        tempos_lqr1=[tempos_lqr1 t1_lqr_span];
        %[t1_lqr,x1_lqr] = ode45(@(t,x)ModuloLunar(t,x,ControleU(:,i)),[T_valido(i) T_valido(i+inc)],y_valido(i,:));
        [t1_lqr,x1_lqr] = ode45(@(t,x)ModuloLunar(t,x,ControleU(:,i)),t1_lqr_span,y_valido(i,:));
    end

    %if(T_valido(i)>=7.37 && T_valido(i)<2880+7.37)
    if(T_valido(i)>=7.37 && T_valido(i)<Tf0+7.37)
        t2_lqr_span = linspace(T_valido(i),T_valido(i+inc),2);
        tempos_lqr1=[tempos_lqr1 t2_lqr_span];
        %[t2_lqr,x2_lqr] = ode45(@(t,x)ModuloLunar(t,x,ControleU(:,i)),[T_valido(i) T_valido(i+inc)],y_valido(i,:));
        [t2_lqr,x2_lqr] = ode45(@(t,x)ModuloLunar(t,x,ControleU(:,i)),t2_lqr_span,y_valido(i,:));
    end
    
    %if(T_valido(i)>=7.37+2880 && T_valido(i)<3100)
    if(T_valido(i)>=7.37+Tf0 && T_valido(i)<Tf1)
        t3_lqr_span = linspace(T_valido(i),T_valido(i+inc),2);
        tempos_lqr1=[tempos_lqr1 t3_lqr_span];
       % [t3_lqr,x3_lqr] = ode45(@(t,x)ModuloLunar(t,x,ControleU(:,i)),[T_valido(i) T_valido(i+inc)],y_valido(i,:));
        [t3_lqr,x3_lqr] = ode45(@(t,x)ModuloLunar(t,x,ControleU(:,i)),t3_lqr_span,y_valido(i,:));
    end

    %if(T_valido(i)>=3100 && T_valido(i)<3400)
    if(T_valido(i)>=Tf1 && T_valido(i)<Tf2)
        t4_lqr_span = linspace(T_valido(i),T_valido(i+inc),2);
        tempos_lqr1=[tempos_lqr1 t4_lqr_span];
        %[t4_lqr,x4_lqr] = ode45(@(t,x)ModuloLunar(t,x,ControleU(:,i)),[T_valido(i) T_valido(i+inc)],y_valido(i,:)); 
        [t4_lqr,x4_lqr] = ode45(@(t,x)ModuloLunar(t,x,ControleU(:,i)),t4_lqr_span,y_valido(i,:)); 
               %[t2_lqr,x2_lqr] = ode45(@(t,x)ModuloLunar(t,x,ControleU(:,i)),[3100 3400],x1_lqr(end,:)); 
    end
    
    if(T_valido(i)>=Tf2 && T_valido(i)<Tf)
        t5_lqr_span = linspace(T_valido(i),T_valido(i+inc),2);
        tempos_lqr1=[tempos_lqr1 t5_lqr_span];
        %[t5_lqr,x5_lqr] = ode45(@(t,x)ModuloLunar(t,x,ControleU(:,i)),[T_valido(i) T_valido(i+inc)],y_valido(i,:)); 
       [t5_lqr,x5_lqr] = ode45(@(t,x)ModuloLunar(t,x,ControleU(:,i)),t5_lqr_span,y_valido(i,:)); 
             %[t5_lqr,x5_lqr] = ode45(@(t,x)ModuloLunar(t,x,ControleU(:,i)),[3400 3510],x2_lqr(end,:)); 
    end
%     if(T_valido(i)>=Tf && T_valido(i)<Tff)
%         [t6_lqr,x6_lqr] = ode45(@(t,x)ModuloLunar(t,x,ControleU(:,i)),[T_valido(i) T_valido(i+inc)],y_valido(i,:)); 
%        %[t5_lqr,x5_lqr] = ode45(@(t,x)ModuloLunar(t,x,ControleU(:,i)),[3400 3510],x2_lqr(end,:)); 
%     end
end

% Reunindo soluções do LQR
t_lqr=[t1_lqr;t2_lqr;t3_lqr;t4_lqr;t5_lqr];%t6_lqr];
x_lqr=[x1_lqr;x2_lqr;x3_lqr;x4_lqr;x5_lqr];%x6_lqr];
% 
% % Cálculo de erros entre o sistema com controlador e a malha aberta
% erro_r_lqr1=[];
% erro_theta_lqr1=[];
% erro_phi_lqr1=[];
% erro_vr_lqr1=[];
% erro_vtheta_lqr1=[];
% erro_vphi_lqr1=[];
% erro_m_lqr1=[];
% 
% for i=1:1:length(T_valido)
%   for j=1:1:length(t_lqr)
%     if(T_valido(i)==t_lqr(j))
%         erro_r_lqr1=[erro_r_lqr1;(y_valido(i,1)-x_lqr(j,1))^2];
% %         erro_theta_lqr1=[erro_theta_lqr1;(y_valido(i,2)-x_lqr(j,2))^2];
% %         erro_phi_lqr1=[erro_phi_lqr1;(y_valido(i,3)-x_lqr(j,3))^2];
% %         erro_vr_lqr1=[erro_vr_lqr1;(y_valido(i,4)-x_lqr(j,4))^2];
% %         erro_vtheta_lqr1=[erro_vtheta_lqr1;(y_valido(i,5)-x_lqr(j,5))^2];
% %         erro_vphi_lqr1=[erro_vphi_lqr1;(y_valido(i,6)-x_lqr(j,6))^2];
% %         erro_m_lqr1=[erro_m_lqr1;(y_valido(i,7)-x_lqr(j,7))^2];     
%     end
%   end
% end


%% Gráficos do sistema controlado vs malha aberta 
figure
plot(t_lqr,x_lqr(:,1),'r-')
hold on
plot(T_PLOT,Y_PLOT(:,1),'b-')
legend('Distância radial com LQR','Distância radial em malha aberta')
xlabel('t(s)')
ylabel('r (km)')
% 
figure
plot(t_lqr,x_lqr(:,2)*180/pi,'r-')
hold on
plot(T_PLOT,Y_PLOT(:,2)*180/pi,'b-')
legend('Longitude com LQR','Longitude em malha aberta')
xlabel('t(s)')
ylabel('theta (°)')
% 
figure
plot(t_lqr,x_lqr(:,3)*180/pi,'r')
hold on
plot(T_PLOT,Y_PLOT(:,3)*180/pi,'b-')
legend('Latitude com LQR','Latitude em malha aberta')
xlabel('t(s)')
ylabel('phi (°)')
%
figure
plot(t_lqr,x_lqr(:,4)*1000,'r')
hold on
plot(T_PLOT,Y_PLOT(:,4)*1000,'b-')
legend('Velocidade radial com LQR','Velocidade radial em malha aberta')
xlabel('t(s)')
ylabel('vr (m/s)')
% 
figure
plot(t_lqr,x_lqr(:,5)*1000,'r')
hold on
plot(T_PLOT,Y_PLOT(:,5)*1000,'b-')
legend('Velocidade longitudinal com LQR','Velocidade longitudinal em malha aberta')
xlabel('t(s)')
ylabel('vtheta (m/s)')
% 
figure
plot(t_lqr,x_lqr(:,6)*1000,'r')
hold on
plot(T_PLOT,Y_PLOT(:,6)*1000,'b-')
legend('Velocidade latitudinal com LQR','Velocidade latitudinal em malha aberta')
xlabel('t(s)')
ylabel('vphi (m/s)')
% 
figure
plot(t_lqr,x_lqr(:,7),'r')
hold on
plot(T_PLOT,Y_PLOT(:,7),'b-')
legend('Fluxo de massa com LQR','Fluxo de massa em malha aberta')
xlabel('t(s)')
ylabel('m (kg)')

%% Função que representa o sistema do Módulo Lunar 
function dxdt = ModuloLunar(t,x,u_controle)
alpha=u_controle(1);
betha=u_controle(2);
k=u_controle(3);

g0=9.81e-3; % gravidade da Terra em Km/s^2
Isp=316; % Impulso Específico em segundos
Tmax=1700e-3; % Força de Propulsão Máxima em Kg*Km/s^2
mi=4902.78; % Parâmetro Gravitacional Padrão em km^3/s^2
omega=2.6632e-6; % Velocidade Angular da Lua em rad/s.


dxdt=[x(4); 
      x(5)/(x(1)*cos(x(3))); 
      x(6)/x(1); 
      ((-Tmax*k)/x(7))*sin(betha) - (mi)/(x(1)^2) + (x(6)^2)/x(1) + (x(5)^2)/x(1) + x(1)*(omega^2)*cos(x(3)) + 2*omega*x(5)*cos(x(3));
      ((Tmax*k)/x(7))*cos(betha)*cos(alpha) - (x(4)*x(5))/x(1) + (x(5)*x(6)*sin(x(3)))/(x(1)*cos(x(3))) + 2*omega*x(6)*sin(x(3)) - 2*omega*x(4)*cos(x(3));
      ((Tmax*k)/x(7))*cos(betha)*sin(alpha) - (x(4)*x(6))/x(1) - ((x(5)^2)*sin(x(3)))/(x(1)*cos(x(3))) - x(1)*(omega^2)*sin(x(3))*cos(x(3)) - 2*omega*x(5)*sin(x(3));
      -(Tmax*k)/(Isp*g0)];

end
