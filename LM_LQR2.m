%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Projeto Final do curso de Engenharia de Controle e Automação
% Universidade: CEFET - RJ/Uned NI
% Aluna: Laís Lima - Matrícula: 1620368ECAN
% Professor orientador: Mauro Vasconcellos
% Referência principal: Artigo "Three-Dimensional Trajectory Optimization of Soft Lunar Landings from the Parking Orbit with Considerations of the Landing Site" escrito por Bong-Gyun Park and Min-Jea Tahk (2011)
% Script: LM_LQR2.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%close all

% Criando variáveis auxiliares
ControleU=[];
y_valido=[R_valido THETA_valido PHI_valido VR_valido VTHETA_valido VPHI_valido M_valido];
y_transposto=y_valido';
estadosX=[];
Matrizes_A=[];
Matrizes_B=[];
u_valido=[ALPHA_valido BETHA_valido K_valido];
u_transposto=u_valido';
entradasU=[];

ref=[1737.4; 0.3546;0.3575;-1e-3;0;0;250]; %estado de referência/estado final

%Reescrevendo  os vetores de estados em uma matriz (7 x m), em que 
%m=length(y_válido). Apenas os estados válidos do sistema foram
%considerados, de acordo com o critério de controlabilidade

for i=1:1:length(y_valido)
    estadosX=[estadosX; y_transposto(:,i)];
end

for i=1:1:length(u_valido)
    entradasU=[entradasU; u_transposto(:,i)];
end

%%% Matrizes Q e R pelo método de Bryson
%tenho que considerar as faixas
Q=zeros(7,7); 
R=zeros(3,3);
%vetor_desvioX=[18 12 4 38 1600 340 350]; %pensar em algum padrão pra ele ir diminuindo os desvios conforme avança rumo à superfície
%vetor_desvioU=[9 30 1];
% 
% j=1;
% for i=1:1:length(vetor_desvioX)
%     Q(i,j)=1/(vetor_desvioX(i))^2;
%     j=j+1;
% end
% j=1;
% for i=1:1:length(vetor_desvioU)
%     R(i,j)=1/(vetor_desvioU(i))^2;
%     j=j+1;
% end

vetor_desvioX=[];
vetor_desvioU=[];

% for i=1:1:length(R_valido)-1
%     atual=R_valido(i);
%     proximo=R_valido(i+1);
%     if proximo>atual
%         desvio_x1=proximo-atual;    
%     else
%         desvio_x1=atual-proximo;
%     end
% end

%Verificando positividade das matrizes Q e R    
% A matriz Q é semi-definida positiva (autovalores_na_diagonal>=0)
positividadeQ=eig(Q);
% A matriz R é definida positiva (autovalores_na_diagonal>0)
positividadeR=eig(R);

% Calculando as leis de controle apenas para a Powered Descent Phase
inc=2;
l=0;
percorreT=1;
for i=1:7:(length(A_valida)-(7*inc+6))
        A_atual=A_valida(i:i+6,:);
        B_atual=B_valida(i:i+6,:);
        ref0=estadosX(i+7*inc:i+7*inc+6,:);
        if(estadosX(i)>ref0(1))
            vetor_desvioX(1)=estadosX(i)-ref0(1);
        else
            vetor_desvioX(1)=ref0(1)-estadosX(i);
        end
        if vetor_desvioX(1)==0
            vetor_desvioX(1)=0.000001;
        end
        if(estadosX(i+1)>ref0(2))
            vetor_desvioX(2)=estadosX(i+1)-ref0(2);
        else
            vetor_desvioX(2)=ref0(2)-estadosX(i+1);
        end
        if vetor_desvioX(2)==0
            vetor_desvioX(2)=0.000001;
        end
        if(estadosX(i+2)>ref0(3))
            vetor_desvioX(3)=estadosX(i+2)-ref0(3);
        else
            vetor_desvioX(3)=ref0(3)-estadosX(i+2);
        end
        if vetor_desvioX(3)==0
            vetor_desvioX(3)=0.000001;
        end
        if(estadosX(i+3)>ref0(4))
            vetor_desvioX(4)=estadosX(i+3)-ref0(4);
        else
            vetor_desvioX(4)=ref0(4)-estadosX(i+3);
        end
        if vetor_desvioX(4)==0
            vetor_desvioX(4)=0.000001;
        end
        if(estadosX(i+4)>ref0(5))
            vetor_desvioX(5)=estadosX(i+4)-ref0(5);
        else
            vetor_desvioX(5)=ref0(5)-estadosX(i+4);
        end
        if vetor_desvioX(5)==0
            vetor_desvioX(5)=0.000001;
        end
        if(estadosX(i+5)>ref0(6))
            vetor_desvioX(6)=estadosX(i+5)-ref0(6);
        else
            vetor_desvioX(6)=ref0(6)-estadosX(i+5);
        end
        if vetor_desvioX(6)==0
            vetor_desvioX(6)=0.000001;
        end
        if(estadosX(i+6)>ref0(7))
            vetor_desvioX(7)=estadosX(i+6)-ref0(7);
        else
            vetor_desvioX(7)=ref0(7)-estadosX(i+6);
        end
        if vetor_desvioX(7)==0
            vetor_desvioX(7)=0.000001;
        end
        %vetor_desvioX=[ref0-estadosX(i) ref0-estadosX(i+1) ref0-estadosX(i+2) ref0-estadosX(i+3) ref0-estadosX(i+4) ref0-estadosX(i+5) ref0-estadosX(i+6)]; %pensar em algum padrão pra ele ir diminuindo os desvios conforme avança rumo à superfície
        
       % vetor_desvioU=[9 30 1];
               
        if(T_valido(percorreT)>=0 && T_valido(percorreT)<7.37)
            k_desvioU=1;
            alpha_desvioU=0.000001;
            betha_desvioU=(-205+(25/7.37)*7.37)*(pi/180); % onde t = 7.37
        end
        if(T_valido(percorreT)>=7.37 && T_valido(percorreT)<Tf0+7.37)
            k_desvioU=0.000001;
            alpha_desvioU=0.000001;
            betha_desvioU=-180*(pi/180);
        end
        if(T_valido(percorreT)>=7.37+Tf0 && T_valido(percorreT)<Tf1)
            k_desvioU=1;
            alpha_desvioU=(6.5*(Tf1-Tf0)/(Tf1-Tf0))*pi/180; % onde t = Tf1
            betha_desvioU=(-180+30*(Tf1-Tf0)/(Tf-Tf0))*pi/180; % onde t = Tf1
        end
        if(T_valido(percorreT)>=Tf1 && T_valido(percorreT)<Tf2)
            k_desvioU=1;
            betha_desvioU=(-180+30*(Tf2-Tf0)/(Tf-Tf0))*pi/180; % onde t = Tf2
            alpha_desvioU=(6.5-0.5*(Tf2-Tf1)/(Tf2-Tf1))*pi/180; % onde t = Tf2
        end
        if(T_valido(percorreT)>=Tf2 && T_valido(percorreT)<Tf)
            k_desvioU=1;
            betha_desvioU=(-180+30*(Tf-Tf0)/(Tf-Tf0))*pi/180; % onde t = Tf
            alpha_desvioU=(6+3*(Tf-Tf2)/(Tf-Tf2))*pi/180; % onde t = Tf
        end
        if(T_valido(percorreT)>=Tf && T_valido(percorreT)<Tff)
            k_desvioU=0.000001;
            alpha_desvioU=0.000001;
            betha_desvioU=0.000001;
        end
        percorreT=percorreT+1;
        vetor_desvioU=[alpha_desvioU betha_desvioU k_desvioU];
%         u0=entradasU(i+3*inc:i+3*inc+2,:);
%         vetor_desvioU=[u0-entradasU(i) u0-entradasU(i+1) u0-entradasU(i+3)];
         z=1;
        for w=1:1:length(vetor_desvioX)
            Q(w,z)=1/(vetor_desvioX(w))^2;
            z=z+1;
        end
        z=1;
        for w=1:1:length(vetor_desvioU)
            R(w,z)=1/(vetor_desvioU(w))^2;
            z=z+1;
        end
        % Calculando o ganho K
        GanhoK_atual=lqr(A_atual,B_atual,Q,R);
        % Montando a lei de controle LQR
        
        ControleU_atual=-GanhoK_atual*(estadosX(i:i+6,:)-ref0);   
        ControleU=[ControleU ControleU_atual];
end

for i=1:inc:(length(ControleU)-(inc+1))

    if(T_valido(i)>=0 && T_valido(i)<7.37)
        [t1_lqr2,x1_lqr2] = ode45(@(t,x)ModuloLunar(t,x,ControleU(:,i)),[T_valido(i) T_valido(i+inc)],y_valido(i,:));
    end

    %if(T_valido(i)>=7.37 && T_valido(i)<2880+7.37)
    if(T_valido(i)>=7.37 && T_valido(i)<Tf0+7.37)
        [t2_lqr2,x2_lqr2] = ode45(@(t,x)ModuloLunar(t,x,ControleU(:,i)),[T_valido(i) T_valido(i+inc)],y_valido(i,:));
    end
    
    %if(T_valido(i)>=7.37+2880 && T_valido(i)<3100)
    if(T_valido(i)>=7.37+Tf0 && T_valido(i)<Tf1)
        [t3_lqr2,x3_lqr2] = ode45(@(t,x)ModuloLunar(t,x,ControleU(:,i)),[T_valido(i) T_valido(i+inc)],y_valido(i,:));
    end

    %if(T_valido(i)>=3100 && T_valido(i)<3400)
    if(T_valido(i)>=Tf1 && T_valido(i)<Tf2)
        [t4_lqr2,x4_lqr2] = ode45(@(t,x)ModuloLunar(t,x,ControleU(:,i)),[T_valido(i) T_valido(i+inc)],y_valido(i,:)); 
        %[t2_lqr,x2_lqr] = ode45(@(t,x)ModuloLunar(t,x,ControleU(:,i)),[3100 3400],x1_lqr(end,:)); 
    end
    
%     if(T_valido(i)>=Tf2 && T_valido(i)<Tf)
%         [t5_lqr,x5_lqr] = ode45(@(t,x)ModuloLunar(t,x,ControleU(:,i)),[T_valido(i) T_valido(i+inc)],y_valido(i,:)); 
%        %[t5_lqr,x5_lqr] = ode45(@(t,x)ModuloLunar(t,x,ControleU(:,i)),[3400 3510],x2_lqr(end,:)); 
%     end
end

% Reunindo soluções do LQR
t_lqr2=[t1_lqr2;t2_lqr2;t3_lqr2;t4_lqr2];%t5_lqr];
x_lqr2=[x1_lqr2;x2_lqr2;x3_lqr2;x4_lqr2];%x5_lqr];

%% Gráficos do sistema controlado vs malha aberta 
figure
plot(t_lqr2,x_lqr2(:,1),'g-')
hold on
plot(T_PLOT,Y_PLOT(:,1),'b-')
legend('Distância radial com LQR','Distância radial em malha aberta')
xlabel('t(s)')
ylabel('r (km)')

figure
plot(t_lqr2,x_lqr2(:,2)*180/pi,'g-')
hold on
plot(T_PLOT,Y_PLOT(:,2)*180/pi,'b-')
legend('Longitude com LQR','Longitude em malha aberta')
xlabel('t(s)')
ylabel('theta (°)')

figure
plot(t_lqr2,x_lqr2(:,3)*180/pi,'g-')
hold on
plot(T_PLOT,Y_PLOT(:,3)*180/pi,'b-')
legend('Latitude com LQR','Latitude em malha aberta')
xlabel('t(s)')
ylabel('phi (°)')

figure
plot(t_lqr2,x_lqr2(:,4)*1000,'g-')
hold on
plot(T_PLOT,Y_PLOT(:,4)*1000,'b-')
legend('Velocidade radial com LQR','Velocidade radial em malha aberta')
xlabel('t(s)')
ylabel('vr (m/s)')

figure
plot(t_lqr2,x_lqr2(:,5)*1000,'g-')
hold on
plot(T_PLOT,Y_PLOT(:,5)*1000,'b-')
legend('Velocidade longitudinal com LQR','Velocidade longitudinal em malha aberta')
xlabel('t(s)')
ylabel('vtheta (m/s)')

figure
plot(t_lqr2,x_lqr2(:,6)*1000,'g-')
hold on
plot(T_PLOT,Y_PLOT(:,6)*1000,'b-')
legend('Velocidade latitudinal com LQR','Velocidade latitudinal em malha aberta')
xlabel('t(s)')
ylabel('vphi (m/s)')

figure
plot(t_lqr2,x_lqr2(:,7),'g-')
hold on
plot(T_PLOT,Y_PLOT(:,7),'b-')
legend('Fluxo de massa com LQR','Fluxo de massa em malha aberta')
xlabel('t(s)')
ylabel('m (kg)')

%%%%%%%%%%%% Função que representa o sistema do Módulo Lunar %%%%%%%%%%%%%

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

