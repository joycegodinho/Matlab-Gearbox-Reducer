clear all
close all
clc

% 
% 
% A = ones(3,5);
% D = A.*3;
% 
% B = [5 6 7 8 9
%      1 2 3 4 5
%      6 5 4 3 2];
% C = D.*B;
m = [4.5 5 6 8 10 12 16 20 25 32 40 50]; %módulos comerciais relevantes
phi = 20;
psi = 30;
Ns = 17:150;
F = (25/24.5)/cosd(30); %Largura de Face.
HB = 334; %Dureza Brinel 9310
A = ones(12,134);

% Fator de sobrecarga K0
        K0 = 1;                       % Choques leves motor elétrico
        
        % Fator Dinamico Kv
        Qv = 10;
        B_ = 0.25*(12-Qv)^(2/3);
        A_ = 50 + 56*(1-B_);
for i = 1:length(m)
    for j = 1:length(Ns)
    H = 600*0.98; %Potencial total do motor

red = 3; %Redução geral do sistema
Na = Ns*(red-1);
Ds(i,j) = Ns(j).*m(i);
Da(i,j) = Na(j).*m(i);
Dp(i,j) = (Da(i,j)-Ds(i,j))./2;
N_p(i,j) = Dp(i,j)./m(i);
N_p = ceil(abs(N_p(1,:)));
red_p(j) = N_p(j)/Ns(j); %Redução para projeto
eng_p = 3; % N# de eng planetárias
He = H/3; %Potencia transmitida em cada engrenamento planeta-solar

n_solar = 13500;
n_planeta = 13500*red_p;
n_braco = 13500/red;

%% Fator de forma de Lewis
    Ng = N_p;
    Np = Ns;
    

if Np<=12                       
        Yp=0.245;
    elseif Np==13
        Yp=0.261;
    elseif Np==14
        Yp=0.277;
    elseif Np==15
        Yp=0.290;
    elseif Np==16
        Yp=0.296;
    elseif Np==17
        Yp=0.303;
    elseif Np==18    
        Yp=0.309;
    elseif Np==19
        Yp=0.314;
    elseif Np==20
        Yp=0.322;
    elseif Np==21
        Yp=0.328;
    elseif Np==22
        Yp=0.331;    
    elseif Np==24
        Yp=0.337;
    elseif Np==26
        Yp=0.346; 
    elseif Np==28
        Yp=0.353;
    elseif Np==30
        Yp=0.359;
    elseif Np==33
        Yp=0.371;
    elseif Np==36
        Yp=0.377;    
    elseif Np==38
        Yp=0.384;
    else
        Yp=0.397;
    end
    
    % Coroa
    if Ng==15
        Yg=0.29;
    elseif Ng==16
        Yg=0.296;
    elseif Ng==17
        Yg=0.303;
    elseif Ng==18    
        Yg=0.309;
    elseif Ng==19
        Yg=0.314;
    elseif Ng==20
        Yg=0.322;
    elseif Ng==21
        Yg=0.328;
    elseif Ng==22
        Yg=0.331;    
    elseif Ng==24
        Yg=0.337;
    elseif Ng==26
        Yg=0.346;
    elseif Ng==27
        Yg=0.353;
    elseif Ng==33
        Yg=0.359;
    elseif Ng==36
        Yg=0.377;
    elseif Ng==39
        Yg=0.384;
    elseif Ng==45
        Yg=0.401;
    elseif Ng==48
        Yg=0.409;
    elseif Ng==60
        Yg=0.422;
    elseif Ng==72
        Yg=0.430;
    elseif Ng==80
        Yg=0.435;
    elseif Ng==84
        Yg=0.438;
    elseif Ng==96
        Yg=0.447;
    elseif Ng==108
        Yg =0.450;
    elseif Ng==120
        Yg =0.454;
    elseif Ng==150
        Yg=0.46;
    else 
        Yg=0.50;
    end
     
     B(i,j) = Np(j).*A(i,j);
     Dsi = Ds./25.4; %metrico para ingles
     Dpi = Dp./25.4;
     P = B./Dsi; % Passos diametrais
     P = P(:,1)'; 
     p = pi*m;
     
      V = pi*Dsi*n_solar/12;                    % Velocidade da linha primitiva, ft/min
      V_m = V.*0.0051; %Velocidade da linha primitiva em m/s;
      Wt = 33000*He./V;                     % Força transversal, em lbf
      W = Wt/(cosd(phi)*cosd(psi));
      Wa = W*cosd(phi)*sind(psi);
      Wr = W*sind(phi);
    % Determinacao dos K
    
        % Fator Dinamico Kv
        Kv(i,j) = ((A_ + (V(i,j)^(1/2)))/A_)^B_;
        
        % Fator espessura de aro Kb
        
        Kb = 1;                          % Razao mb>1.2
    
        % Fator de Confiabilidade Kr
        
        conf = 0.99;                      % Confiabilidade
        Kr = 0.5 - 0.109*log(1-conf);
          
        % Fator Ciclagem de Tensao Yn e Zn para 10^7 ciclos
        
            % Fatores para o Pinhao
            Ynp = 1;            
            Znp = 1;
            
            % Fatores para a coroa
            Yng = 1;              
            Zng = 1;
            
        % Fator Tamanho Ks
        
            Ksp = 1.192.*((F.*sqrt(Yp))./P).^0.0535; % Para o pinhao
            Ksg = 1.192.*((F.*sqrt(Yg))./P).^0.0535; % Para a coroa
            
        % Fator Distribuicao de Carga Km
        
            Cmc = 1;            % Dentes sem coroamento
            Cpm = 1.1;          % Valor mais conservador, nao necessariamente as engrenagens estarao centralizadas entre mancais
            Ce = 1;             % Outras condiçoes
            A_m = 0.0675;       % Engrenamento preciso
            B_m = 0.0128;       % Engrenamento preciso
            C_m = -0.0926e-4;   % Engrenamento preciso
            
            Cma = A_m +(B_m*F)+(C_m*F^2);
            
            if (F<=1)
            Cpf_p = (F./(10.*Dsi))-0.025;
            Cpf_g = (F./(10.*Dpi))-0.025;
            
            elseif ((F>1)&&(F<=17))
            Cpf_p = (F./(10.*Dsi))-0.0375+(0.0125*F);
            Cpf_g = (F./(10.*Dpi))-0.0375+(0.0125*F);
            
            else
            Cpf_p = (F./(10.*Dsi))-0.1109+(0.0207*F)-(0.000228*F^2);
            Cpf_g = (F./(10.*Dpi))-0.1109+(0.0207*F)-(0.000228*F^2);
            end
            
            Kmp = 1 + Cmc*(Cpf_p*Cpm + Cma*Ce); % Para o pinhao
            Kmg = 1 + Cmc*(Cpf_g*Cpm + Cma*Ce); % Para a coroa
            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % COEFICIENTE ELÁSTICO CP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TABELA 14.8
    Cp = 2300; %sqrt(psi)
    
    % Determinacao do J, Fator Geometrico de Resistencia a Flexao
    
        if Np==14
        Jp = 0.23;
        Jg = 0.37;
        elseif Np==15
        Jp = 0.25;
        Jg = 0.38;
        elseif Np==16
        Jp = 0.27;
        Jg = 0.39;
        else  Jp=0.33;
        Jg = 0.25;
        end
        
    % Determinacao do I, Fator Geometrico de Resistencia Superficial
    
        mn = 1;                                         % Engrenagens cilindricas de dentes retos
        Ip =(cosd(phi)*sind(phi).*red_p)./(2*mn.*(red_p + 1));   % Para pinhao, engrenagem externa
        Ig =(cosd(phi)*sind(phi).*red_p)./(2*mn.*(red_p + 1));   % Para coroa, engrenagem externa
        
    % Fator de temperatura Kt
        Kt = 1;
        
    % Fator de condicao de superfice Cf
        Cf = 1;
       
    % Fator razao de Dureza Ch
    
        %if hbsp/hbsg > 1.7
         %   A1 = 0.00698;
        %elseif hbsp/hbsg < 1.2
         %   A1 = 0;
        %else 
        %    A1 = 8.98*10^-3*(hbsp/hbsg)-8.29*10^-3;
        %end
        Ch = 1;      


    % Calculo das tensoes
    
    % a) Tensao de flexao St
    % AÇO 9310 
    Stp = 102*HB+16400; %grafico pagina 730-731, figura 14.4 %Inserir em psi
    Stg = 102*HB+16400;
  
    
    % b) Tensao de contato Sc
    % ACO 9310
    Scp = 349*HB+34300;%Resistência ao contato - tabela 14-6 %inserir em psi
    Scg = 349*HB+34300;
    Ksp(i,j) = Ksp(i).*A(i,j);
    Ksg(i,j) = Ksg(i).*A(i,j);
    %B(i,j) = Np(j).*A(i,j);
     P = B./Dsi;
     Ip(i,j) = Ip(j).*A(i,j);
     Ig(i,j) = Ig(j).*A(i,j);
    %Tensão de Flexão da Engrenagem [psi]
    Sigma_fp(i,j) = Wt(i,j).*K0.*Kv(i,j).*Ksp(i,j).*P(i,j).*Kmp(i,j).*Kb./(F*Jp); 
    Sigma_fg(i,j) = Wt(i,j).*K0.*Kv(i,j).*Ksg(i,j).*P(i,j).*Kmg(i,j).*Kb./(F*Jg);

    % Fator de seguranca para Flexao
    Sfp(i,j) = Stp.*Ynp./(Kt.*Kr.*Sigma_fp(i,j));
    Sfg(i,j) = Stg.*Zng./(Kt.*Kr.*Sigma_fg(i,j));

    % Tensao de contato da Engrenagem
    Sigma_cp(i,j) = Cp.*sqrt(Wt(i,j).*K0.*Kv(i,j).*Ksp(i,j).*Kmp(i,j).*Cf./(Dsi(i,j).*F.*Ip(i,j))); % Eq 14-16
    Sigma_cg(i,j) = Cp.*sqrt(Wt(i,j).*K0.*Kv(i,j).*Ksg(i,j).*Kmg(i,j).*Cf./(Dpi(i,j).*F.*Ig(i,j)));
    
    % Fator de seguranca para desgaste
    Shp(i,j) = Scp.*Znp./(Kt.*Kr.*Sigma_cp(i,j)); % Eq 14-42
    Shg(i,j) = Scg.*Zng.*Ch./(Kt.*Kr.*Sigma_cg(i,j));
    
    
    end 
end

[col] = find (Shg(1,:)>1 & Shg(1,:)<1.3, 100);
for k = 1:length(col)
    coluna = col(k);
    SHG = Shg(1,coluna);
    WT = Wt(1,coluna);
    WA = Wa(1,coluna);
    WR = Wr(1,coluna);
    N_s = Ns(coluna);
    NP = N_p(coluna);
    NA = Na(coluna);
    if SHG > 1.15;
        Fator_contato = SHG;
        Dentes_Solar = N_s;
        Dentes_Planeta = NP;
        Dentes_Anelar = NA;
        Diametro_Solar = Dentes_Solar*4.5;
        Diametro_Anelar = Dentes_Anelar*4.5;
        Diametro_Planeta = Dentes_Planeta*4.5;
        Forca_Tangencial = WT*4.44822;
        Forca_Axial = WA*4.44822;
        Forca_Radial = WR*4.44822;
    end
end
    
disp(['Fator_contato = ' num2str(Fator_contato)])
disp(['Dentes_Solar = ' num2str(Dentes_Solar)])
disp(['Dentes_Anelar = ' num2str(Dentes_Anelar)])
disp(['Dentes_Planeta = ' num2str(Dentes_Planeta)])
disp(['modulo = 4.5 '])
disp(['Diametro_Solar = ' num2str(Diametro_Solar) '[mm]'])
disp(['Diametro_Anelar = ' num2str( Diametro_Anelar) '[mm]'])
disp(['Diametro_Planeta = ' num2str(Diametro_Planeta) '[mm]'])
disp(['Forca_Tangencial = ' num2str(Forca_Tangencial) '[N]'])
disp(['Forca_Axial = ' num2str(Forca_Axial) '[N]'])
disp(['Forca_Radial = ' num2str(Forca_Radial) '[N]'])




% Red_P = find(red_p(col)>0.4 & red_p(col)<0.6, 12);
% RedP = red_p(Red_P);
% N_Solar = find(Ns(col)>17 & Ns(col)<150, 12);
% Dentes_Solar = Ns(N_Solar);
% M = find(m(row)>4.5 & m(row)<20, 12);
% Modulo = m(M);
% Dentes_Anelar = (red-1)*Dentes_Solar;
% Dentes_Planeta = floor(Dentes_Solar.*(RedP));