%%Engrenamento Simples
clc; clear all; close all

%%TEP 1
%%Dados de entrada
n_entrada = 13500; %rpm
n_saida = n_entrada/3; %rpm
H_entrada = 600; %hp
H_entrada = H_entrada/3; %Devido ao número de engrenagens planetárias
reducao = 3;

%%Fatores de Projeto
nd = 2; %fator de projeto
Pd = 12; %passo Diametral
R = 0.99; %confiabilidade
k1 = 1; %dentes completos
%%Definir Número de Dentes do TEP

%Assumindo
phi_normal = 20 *pi/180;
phi_helice = 30 *pi/180;
phi_transversal = (atan(tan(phi_normal)/(cos(phi_helice))))*pi/180;
Fator_de_Qualidade = 10;
Qv = Fator_de_Qualidade; 
Lf = 1.5; %Largura de face
%%Definições de Material

Hb = 334;
St = 102*Hb+16400; %grafico pagina 730-731, figura 14.4 %Inserir em psi
Sc = 349*Hb+34300; %Resistência ao contato - tabela 14-6 %inserir em psi

%%Definição de vida -- CHECAR
ciclos = 10^7;
%%Foi selecionado o ciclo de vida mais comum

%%Determinar número de dentes

%%No TEP
%%Condições assumidas
%%1. entrada pela engrenagem central; 2. saída pelo braço
%%3. Solar fixa
%%atencao aos sentidos de rotacao

n_solar = n_entrada; %%solar e braço giram no mesmo sentido
n_braco = n_entrada/reducao;
n_anelar = 0;

%Determinar trem de engrenagem
e = -(n_anelar - n_braco)/(n_solar-n_braco);

N_planetaria = [20:35];
N_central = [13];

        for jj=1:length(N_planetaria)-1
            
N_anelar(jj) = 2*N_planetaria(jj) + N_central;
razao(jj) = -N_central/N_planetaria(jj);
n_planetaria(jj) = razao(jj)*(n_solar-n_braco)+(n_braco);

 %%A tem que ser inteiro para melhor desempenho do TEP
 
A(jj) = (N_central+N_anelar(jj))/3;
        jj= jj+1;
       
        end
       
       nn = [3,6,9,12,15];
        N_central = [18,18,18,18,18];
        N_anelar = N_anelar(nn);
        N_planetaria = N_planetaria(nn);

%%Fator de Lewis - Depende do número de dentes Y(1) = 16 dentes
arquivo = fopen('Fator_Lewis.txt');
Y = fscanf(arquivo,'%f');
fclose(arquivo);


%%Parte1 - Engrenamento entre Planetária e solar
for ii = 1:length (N_central)
    
    %%Definindo Coroa e Pinhão
        dp_planetaria(ii) = N_planetaria(ii)/Pd; %Inc
        dp_central(ii) = N_central(ii)/Pd; %Inc
       
        if dp_planetaria(ii)>dp_central(ii)   
            Ng(ii) = N_planetaria(ii);
            Np(ii) = N_central(ii);
        else
            Np(ii) = N_planetaria(ii);
            Ng(ii) = N_central(ii);
            
        end
            dp_pinhao(ii) = Np(ii)/Pd;
            phi_transversal = atan(tan(phi_normal)/cos(phi_helice));
            Pd = Pd*cos(phi_transversal);

dp_pinhao(ii) = Np(ii)/Pd;
dp_coroa(ii) = Ng(ii)/Pd;

Mg(ii) = Ng(ii)/Np(ii);  %RazÃ£o de Engrenamento
          
Yp = Y(Np(ii) - 11);

%%Fatores Geométricos
 arquivo = fopen('J.txt');
 J = fscanf(arquivo,'%f');
 fclose(arquivo);
 Jp(ii) = J(ii);
 Jg(ii)= 0.62;

 %%Fatores Geométricos para engrenagens de dentes retos
% Jp(ii) = 0.32;
% Jg(ii) = 0.37;

Mn = 1; %RazÃ£o de Compartilhamento de carga 
I_geometrico = ((cos(phi_transversal)*sin(phi_transversal))/2*Mn)*(Mg(ii)/(Mg(ii)+1));%%Engrenagens externas
%%I_geometrico =((cosd(phi_transversal)*sind(phi_transversal))/2*Mn)*(Mg(ii)/(Mg(ii)+1))

%%Lembrar de mudar se mudar a relação entre os dentes
% Yn_p = 1.3558*ciclos^(-0.0178); %input('Figura 14-14'); --MUDAR DEPOIS
% Yn_g = 1.3558*(ciclos/4)^(-0.0178);
% Zn_p = 1.4488*ciclos^(-0.023);
% Zn_g = 1.4488*(ciclos/Pd)^(-0.023);
%%Os fatores de Ciclagem são iguais a 1, pois consederou-se 10^7 Ciclos
Yn_p = 1;
Yn_g = 1;
Zn_p = 1;
Zn_g = 1;

% %%Engrenaagem de dentes retos
% V(ii) = (pi*dp_pinhao(ii)*n_entrada/12); %ft/min %Velocidade de linha primitiva
% Wt(ii) = 33000*H_entrada/V(ii); %lbf; %Equacao 13-35, Shigley 10 Edicao

%%Engrenagem de dentes helicoidais
V(ii) = (pi*dp_pinhao(ii)*n_entrada/12); %ft/min %Velocidade de linha primitiva
Wt(ii) = 33000*H_entrada/V(ii); %lbf; %Equacao 13-35, Shigley 10 Edicao
W(ii) = Wt(ii)/(cos(phi_normal)*cos(phi_helice));
Wr(ii) = W(ii)*sin(phi_normal);
Wa(ii) = W(ii)*cos(phi_normal)*sin(phi_helice);
    
%%Coeficiente elastico
%%Pinhao e coroa de aço Com E = 30.10^6
Cp = 2300; %sqrt(psi)

%%DEFININDO FATORES;
%%Fator Dinâmico
B = 0.25*(12-Qv)^(2/3);
A = 50+56*(1-B);
Kv = ((A+(V(ii)^(1/2)))/A)^B; %Velocidade em ft/min

%Fator de Condição de superfície
Cf = 1;

%%Fator de Tamanho Ks
Ks = 1.192*((Lf*sqrt(Yp))/Pd)^(0.0535);
    if Ks<1;
Ks = 1;
    else
Ks=Ks;
    end
%FATORES A MUDAR
Kb = 1;
Kt = 1;

Cmc = 1; %para dentes sem coroamento
if Lf<1
    if (Lf/(10*dp_pinhao(ii))) <0.05
      Cpf = 0.05-0.025;  
    else   
    Cpf = (Lf/(10*dp_pinhao(ii)))-0.025;
    end
else
Cpf = (Lf/(10^dp_pinhao(ii)))-0.0375+0.0125*Lf; %Supondo Lf<17in
end
Cpm = 1; %assumindo com base no intervalo entre mancais
Cma = (0.0675+0.0128*Lf+((-1)*0.926*(10^(-4))*Lf^2)); %Fechada, de precisÃ£o
Ce = 1; 
%Fator de Distribuição de Carga
Km = 1+Cmc*(Cpf*Cpm+Cma*Ce);
%confiabilidade 
%Kr = 0.658 - 0.0759*log(1-R);
Kr = 1; %Tabela 14-10


%FLEXAO DE ENGRENAGENS DE DENTES RETOS
F_flexao(ii) = nd*Wt(ii)*Kv*Ks*Pd*(Km*Kb/Jp(ii))*(Kr/St*Yn_p);
F_desgaste(ii) = ((Cp*Kr/Sc*Zn_p)^2)*(nd*Wt(ii)*Ks*Kv*Km*Cf/(dp_pinhao(ii)*I_geometrico));

sigma_flexao(ii) = nd*Wt(ii)*Kv*Ks*(Pd/Lf)*(Km*Kb/Jp(ii));
Sflexao_p(ii) = St*Yn_p/sigma_flexao(ii);

sigma_flexao(ii) = nd*Wt(ii)*Kv*Ks*(Pd/Lf)*(Km*Kb/Jg(ii));
Sflexao_g(ii) = St*Yn_p/sigma_flexao(ii);

%%DESGASTE DE ENGRENAGENS DE DENTES RETOS
sigma_desgaste_p(ii) = Cp*(Wt(ii)*Kv*Ks*(Km*Cf/(dp_pinhao(ii)*Lf*I_geometrico)))^(1/2);
Sdesgaste_p(ii) = ((Sc*Zn_p/(Kt*Kr))/sigma_desgaste_p(ii))^2;

sigma_desgaste_g(ii) = Cp*(Wt(ii)*Kv*Ks*(Km*Cf/(dp_coroa(ii)*Lf*I_geometrico)))^(1/2);
Sdesgaste_g(ii) = ((Sc*Zn_g/(Kt*Kr))/sigma_desgaste_g(ii))^2;

ii = ii+1;

end

  ll = [4,5]


    Sdesgaste_p1 = Sdesgaste_p  (ll)
    Sdesgaste_g1 = Sdesgaste_g (ll)
    Sflexao_p1 = Sflexao_p(ll)
    Sflexao_g1 = Sflexao_g(ll)

%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%


%%
%%Parte2 - Engrenamento entre Planetária e Anelar

ii = 0;

for ii = 1:length (N_central)
   
    %%Definindo Coroa e Pinhão
        dp_planetaria(ii) = N_planetaria(ii)/Pd; %Inc
        dp_anelar(ii) = N_anelar(ii)/Pd; %Inc
       
        if dp_planetaria(ii)>dp_anelar(ii)   
            Ng(ii) = N_planetaria(ii);
            Np(ii) = N_anelar(ii);
        else
            Np(ii) = N_planetaria(ii);
            Ng(ii) = N_anelar(ii);
            
        end
            dp_pinhao(ii) = Np(ii)/Pd;
            phi_transversal = atan(tan(phi_normal)/cos(phi_helice));
            Pd = Pd*cos(phi_transversal);

dp_pinhao(ii) = Np(ii)/Pd;
dp_coroa(ii) = Ng(ii)/Pd;

Mg(ii) = Ng(ii)/Np(ii);  %RazÃ£o de Engrenamento
          
Yp = Y(Np(ii) - 11);

%%Fatores Geométricos
 arquivo = fopen('J.txt');
 J = fscanf(arquivo,'%f');
 fclose(arquivo);
 Jp(ii) = J(ii);
 Jg(ii)= 0.62;

 %%Fatores Geométricos para engrenagens de dentes retos
% Jp(ii) = 0.32;
% Jg(ii) = 0.37;

Mn = 1; %RazÃ£o de Compartilhamento de carga 
I_geometrico = ((cos(phi_transversal)*sin(phi_transversal))/2*Mn)*(Mg(ii)/(Mg(ii)+1));%%Engrenagens externas
%%I_geometrico =((cosd(phi_transversal)*sind(phi_transversal))/2*Mn)*(Mg(ii)/(Mg(ii)+1))

%%Lembrar de mudar se mudar a relação entre os dentes
% Yn_p = 1.3558*ciclos^(-0.0178); %input('Figura 14-14'); --MUDAR DEPOIS
% Yn_g = 1.3558*(ciclos/4)^(-0.0178);
% Zn_p = 1.4488*ciclos^(-0.023);
% Zn_g = 1.4488*(ciclos/Pd)^(-0.023);
%%Os fatores de Ciclagem são iguais a 1, pois consederou-se 10^7 Ciclos
Yn_p = 1;
Yn_g = 1;
Zn_p = 1;
Zn_g = 1;

% %%Engrenaagem de dentes retos
% V(ii) = (pi*dp_pinhao(ii)*n_entrada/12); %ft/min %Velocidade de linha primitiva
% Wt(ii) = 33000*H_entrada/V(ii); %lbf; %Equacao 13-35, Shigley 10 Edicao

%%Engrenagem de dentes helicoidais
V(ii) = (pi*dp_pinhao(ii)*n_entrada/12); %ft/min %Velocidade de linha primitiva
Wt(ii) = 33000*H_entrada/V(ii); %lbf; %Equacao 13-35, Shigley 10 Edicao
W(ii) = Wt(ii)/(cos(phi_normal)*cos(phi_helice));
Wr(ii) = W(ii)*sin(phi_normal);
Wa(ii) = W(ii)*cos(phi_normal)*sin(phi_helice);
    
%%Coeficiente elastico
%%Pinhao e coroa de aço Com E = 30.10^6
Cp = 2300; %sqrt(psi)

%%DEFININDO FATORES;
%%Fator Dinâmico
B = 0.25*(12-Qv)^(2/3);
A = 50+56*(1-B);
Kv = ((A+(V(ii)^(1/2)))/A)^B; %Velocidade em ft/min

%Fator de Condição de superfície
Cf = 1;

%%Fator de Tamanho Ks
Ks = 1.192*((Lf*sqrt(Yp))/Pd)^(0.0535);
    if Ks<1;
Ks = 1;
    else
Ks=Ks;
    end
%FATORES A MUDAR
Kb = 1;
Kt = 1;

Cmc = 1; %para dentes sem coroamento
if Lf<1
    if (Lf/(10*dp_pinhao(ii))) <0.05
      Cpf = 0.05-0.025;  
    else   
    Cpf = (Lf/(10*dp_pinhao(ii)))-0.025;
    end
else
Cpf = (Lf/(10^dp_pinhao(ii)))-0.0375+0.0125*Lf; %Supondo Lf<17in
end
Cpm = 1; %assumindo com base no intervalo entre mancais
Cma = (0.0675+0.0128*Lf+((-1)*0.926*(10^(-4))*Lf^2)); %Fechada, de precisÃ£o
Ce = 1; 
%Fator de Distribuição de Carga
Km = 1+Cmc*(Cpf*Cpm+Cma*Ce);
%confiabilidade 
%Kr = 0.658 - 0.0759*log(1-R);
Kr = 1; %Tabela 14-10


%FLEXAO DE ENGRENAGENS DE DENTES RETOS
sigma_flexao(ii) = nd*Wt(ii)*Kv*Ks*(Pd/Lf)*(Km*Kb/Jp(ii));
Sflexao_p(ii) = St*Yn_p/sigma_flexao(ii);

sigma_flexao(ii) = nd*Wt(ii)*Kv*Ks*(Pd/Lf)*(Km*Kb/Jg(ii));
Sflexao_g(ii) = St*Yn_p/sigma_flexao(ii);

%%DESGASTE DE ENGRENAGENS DE DENTES RETOS
sigma_desgaste_p(ii) = Cp*(Wt(ii)*Kv*Ks*(Km*Cf/(dp_pinhao(ii)*Lf*I_geometrico)))^(1/2);
Sdesgaste_p(ii) = ((Sc*Zn_p/(Kt*Kr))/sigma_desgaste_p(ii))^2;

sigma_desgaste_g(ii) = Cp*(Wt(ii)*Kv*Ks*(Km*Cf/(dp_coroa(ii)*Lf*I_geometrico)))^(1/2);
Sdesgaste_g(ii) = ((Sc*Zn_g/(Kt*Kr))/sigma_desgaste_g(ii))^2;

ii = ii+1;

end


 ll = [4]
 
    N_planetaria = N_planetaria(ll);
    N_anelar =  N_anelar (ll);
    N_central = N_central (ll);
    
    Lf = Lf * 25.4
    
    dp_anelar = dp_anelar(ll) *25.4 %mm
    dp_planetaria = dp_planetaria(ll) *25.4
    dp_central = dp_central(ll) *25.4
    
    densidade_material = 8.0*(10^-6); %kg/mm^3
   
    
%     peso_planetaria = 2*pi*((dp_planetaria/2)^2)*espessura*densidade_material
%       peso_solar = 2*pi*((dp_solar/2)^2)*espessura*densidade_material
%         peso_anelar = 2*pi*((dp_anelar/2)^2)*espessura*densidade_material
%     
    
  
     Sdesgaste_p2 = Sdesgaste_p  (ll)
    Sdesgaste_g2 = Sdesgaste_g (ll)
    Sflexao_p2 = Sflexao_p(ll)
    Sflexao_g2 = Sflexao_g(ll)
    
