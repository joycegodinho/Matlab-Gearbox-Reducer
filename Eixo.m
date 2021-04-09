clear all
close all 
clc
 
% Fator de segurança global para projeto de elementos
n = 1.4; 
 
%%
% DADOS CONHECIDOS
 
% A) EFICIENCIAS DOS MANCAIS E DOS TEPS
 
eta_mancal = 0.98; % Eficiência do par de mancais
eta_tep = 0.95;   % Eficiência do TEP
 
% B) EFICIENCIAS DOS EIXOS
%       Os eixos são orientados da seguinte forma: EIXO 1 é o eixo que
%       conecta o motor a solar do TEP1; EIXO 2 é o eixo que conecta o braço 
%       do TEP1 a solar do TEP2; EIXO 3 é o eixo que conecta o braço do TEP2 a helice.
  
 
 
eta_eixo1 = 1;
eta_eixo2 = eta_mancal*eta_tep;
eta_eixo3 = eta_mancal*eta_tep*eta_mancal*eta_tep;
 
% C) POTÊNCIA NOS EIXOS
 
H_0 = 4.4741994*10^5; 
H_1 = H_0*eta_eixo1;
H_2 = H_1*eta_eixo2;
H_3 = H_2*eta_eixo3;
 
% D) VELOCIDADES ANGULARES E TORQUES NOS EIXOS
%   d.1) Velocidades expressas em rpm. Note que a notação obedece a
%   convenção para os eixos estabelecida em (B).
    n_1 = 13500;
    n_2 = 4500;
    n_3 = 900;
    
%   d.2) Velocidades em rad/s
    omega_1 = (2*pi/60)*n_1;
    omega_2 = (2*pi/60)*n_2;
    omega_3 = (2*pi/60)*n_3;
 
%   d.3) Torques nos eixos, em Nmm
    T_0 =  (H_0/omega_1)*10^3;
    T_1 =  (H_1/omega_1)*10^3;
    T_2 =  (H_2/omega_2)*10^3;
    T_3 =  (H_3/omega_3)*10^3;
 
% E) DADOS RELATIVOS AOS TEPS
    
    g = 9.81;            % Gravidade, m/s^2
    F = (25)/cosd(30); %Largura de Face.
    % e.1) 1o TEP
    Dsol_1 = 157.5 ;       % Diametro da solar, mm
    Dpla_1 = 81;        % Diametro da planeta, mm
    f_1 = 25;           % Largura de face, mm
    Wt_1 = 1312.8256;  % Forca tangencial, N
    Wr_1 = 551.7499;   % Forca radia, N
    Wa_1 = 757.9602;   % forca axial
    Vsol_1 = 3.14*((Dsol_1/2)^2)*F; % Volume solar
    Vpla_1 = 3.14*((Dpla_1/2)^2)*F; % Volume planeta
    Psol_1 = 7861*Vsol_1*10^-9*g;      % Peso da solar, N
    Ppla_1 = 7861*Vpla_1*10^-9*g;     % Peso da planeta, N
    
    
    
    % e.2) 2o TEP
    
    Dsol_2 = 68.9;       % Diametro da solar, mm
    Dpla_2 = 188;        % Diametro da planeta, mm
    f_2 = 176.2;           % Largura de face, mm
    Wt_2 = 11400.224;   % Forca tangencial, N
    Wr_2 = 4790.496;   % Forca radia, N
    Vsol_2 = 3.14*((Dsol_2/2)^2)*F; % Volume solar
    Vpla_2 = 3.14*((Dpla_2/2)^2)*F; % Volume planeta
    Psol_2 = 7861*Vsol_2*10^-9*g;      % Peso da solar, N
    Ppla_2 = 7861*Vpla_2*10^-9*g;     % Peso da planeta, N
    dmancal_1 = 15          % Distância do ultimo mancal do eixo 1 ao TEP1 
    
% F) Considerações de material: AÇO 4340 TROCAR MATERIAL
 
        S_ut = 966;         % Resistencia a ruptura, MPa
        S_y = 897;          % Resistencia ao escoamento, MPa
        S_e1 = 0.5*S_ut;    % Resistencia a fadiga, MPa
        
        % Fatores

        %Fator de Superfície Ka (função da qualidade de acabamento da superfície)
        %k_a = a*S_ut^b  definir acabamento
        k_a = (4.51)*S_ut^(-0.265);  %Acabamento usinado

        %Fator de tamanho Kb
        %Estimar o tamanho do eixo para corrigir esse fator
        k_b = 0.6;

        %Fator de carregamento Kc 
        %corrige o limite de resistência à fadiga quando os ensaios são realizados 
        %com flexão rotativa (1,0), carregamento axial (0,85) e carregamento torcional (0,59). 
        %Como os eixos estão sujeitos a carregamentos combinados, arbitra-se o valor de 0.7.
        k_c = 0.7;

        %Fator de Temperatura Kd
        %Arbitra-se que as temperaturas estarão abaixo de 250 ºC e adota-se Kd = 1. 
        k_d = 1;
        
        %Fator de Confiabilidade Ke
        k_e = 0.814; %Confiabilidade de 99% pois trata-se de um redutor de avião

        %Fator de efeitos diversos Kf
        k_f = 1;
        
        % Resistencia a fadiga, vida infinita.
        S_e = k_a*k_b*k_c*k_d*k_e*k_d*k_f*S_e1;
      
%%
%       1 - DIMENSIONAMENTO DO EIXO DE SAIDA OU EIXO 1
                   clear n
                   n = 2.2
                   Rv = Psol_2; % peso da solar 2
                   Rh = Wt_2*(1-2*cosd(30));
                
                   % DIMENSIONAMENTO NA SEÇÃO CRITICA
                
                   % Na seção critica, não há concetradores de tensão, ou
                   % seja:
                    Kf = 1; % Fator de concentração de tensão para flexao
                    Kfs = 1; % Fator de concentração de tensao para cisalhamento
                  
                  
                    M_a = Ppla_1*dmancal_1;  
                    T_a = 0;
                    M_m = 0
                    T_m = T_1;
                    
                    % Von Mises
                    
                    d_vm = ((16*n/pi)*(1/S_y)*sqrt(4*M_m^2 + 3*T_m^2))^(1/3);
 
                    % Goodman
 
                    A2 = sqrt((4*(Kf*M_a)^2 + 3*(Kfs*T_a)^2));
                    B2 = sqrt((4*(Kf*M_m)^2 + 3*(Kfs*T_m)^2));
 
                    d_Goodman = ((16*n/pi)*((1/S_e )*A2 + (1/S_ut )*B2))^(1/3);
 
                    % Asme Eliptico
 
                    d_ASME = ((16*n/pi)*sqrt( 4*(Kf*M_a/S_e)^2 + ...
       3*(Kfs*T_a/S_e)^2 + 4*(Kf*M_m/S_y)^2 + 3*(Kfs*T_m/S_y)^2))^(1/3);
   
          % DIMENSIONAMENTO POR DEFLEXAO
%       Comprimento do eixo
%       No Eixo 1 são instalados a engrenagem solar do 1o TEP.
 %      
 
               %  Assume-se que o angulo maximo permitido por deflexão para um par de engrenagens é de 0.0005
               % rad. 
              y_all_e = dmancal_1*sin(0.0005);
               
%           
%              Deflexao na direcao y
               
               % A esquerda
               Ay_e =  dmancal_1; % mm
               E = 200*10^3; % MPa
               L_e = dmancal_1;
               d_deflexao_y_e = (abs((n*64*Psol_1*Ay_e^2)*(Ay_e - 3*L_e)/(6*E*pi*y_all_e)))^(1/4);
 

 