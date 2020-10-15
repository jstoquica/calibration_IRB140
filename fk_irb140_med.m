function[Mx My Mz M] = fk_irb140_med(Theta)

%% Revisão
% Em Wang et al. (2018) propoem erros nas juntas (rot) de 0.01 graus
% Em Hafezipour e Khodaygan (2017) propoem erros em cada junta: 0.3-0.4mm 0.2-0.3 graus
% Em Kato et al. (2017) propoem erros em cada junta: 0.3mm (segundo ISO2768) e 0.051 graus
% Taek (2011)consegue medir erros nas juntas de 0.16 graus
% Lind (2012) propoe erros nas juntas de 0.001-0.002 rad (0.057-0.114 graus)


%% Gerar erros sistematicos de comprimento e angulo
% rng(1994,'twister'); % random controller
% eS_c = 2*0.15*(randi([0 1],1,21))-0.15;%segundo ISO2768 (Kato et al., 2017) 0.001% do maior comprimento (+- 0.18mm)
% rng(1994,'twister'); % random controller
% eS_a = 2*0.05*(randi([0 1],1,21))-0.05;% (+- 1.2graus)

% rng(1994,'twister'); % random controller
% eS_c = 0.15*(ones(1,21));%segundo ISO2768 (Kato et al., 2017) 0.001% do maior comprimento (+- 0.18mm)
% rng(1994,'twister'); % random controller
% eS_a = 0.05*(ones(1,21));% (+- 1.2graus)
% 
rng(101994,'twister'); % random controller
eS_c = normrnd(0,0.1,[1,21]);%segundo ISO2768 (Kato et al., 2017) 0.001% do maior comprimento (+- 0.18mm)
rng(111960,'twister'); % random controller
eS_a = normrnd(0,0.2,[1,21]);% (+- 1.2graus)


%% Gerar erros aleaotrios de comprimento e angulo
rng(51989,'twister'); % random controller
eP_c = normrnd(0.1,0.01,[1,18]); % random errors values with normal distribution lenght 
% eP_c = eP_c .* p_sign;
rng(51989,'twister'); % random controller
eP_a = normrnd(0.15,0.01,[1,21]); % random errors values with normal distribution angle
% eP_a = eP_a .* p_sign;


%Parâmetros base
Th0 = 0 + eS_a(1);

%Parâmetros Junta_1
Th1 = Theta(:,1);

%parâmetros Junta_2 - Hayati
% d2 = 0;
Th2 = Theta(:,2);

%Parâmetros Junta_3
Th3 = Theta(:,3);

%Parâmetros Junta_4
Th4 = Theta(:,4);

%Parâmetros Junta_5
Th5 = Theta(:,5);

%Parâmetros Junta_6
Th6 = Theta(:,6);

c1 = cosd(Th1 + eS_a(2));%0.5);
s1 = sind(Th1 + eS_a(2));%0.5);
c2 = cosd(Th2 + eS_a(3));%0.3);
s2 = sind(Th2 + eS_a(3));%0.3);
c3 = cosd(Th3 + eS_a(4));%0.2);
s3 = sind(Th3 + eS_a(4));%0.2);
c4 = cosd(Th4 + eS_a(5));%0.4);
s4 = sind(Th4 + eS_a(5));%0.4);
c5 = cosd(Th5 + eS_a(6));%0.15);
s5 = sind(Th5 + eS_a(6));%0.15);
c6 = cosd(Th6 + eS_a(7));%0.2);
s6 = sind(Th6 + eS_a(7));%0.2);

Kc = 1.15;
Ka = 1.05;

%Parâmetro Base-mundo
alfa0 = 0 + eS_a(8); beta0 = 0 + eS_a(9);
a0 = 0 + eS_c(1); b0 = 0 + eS_c(2); d0 = 0 + eS_c(3);

d1 = 352 + eS_c(4);%0.7;
a1 = 70 + eS_c(5);%0.5;
alfa1 = -90 + eS_a(10);%0.17;
beta1 = 0 + eS_a(11)*Ka;%0.2; %new for full model
b1 = 0 + eS_c(6)*Kc;%0.1; %new for full model
 
d2 = 0 + eS_c(7)*Kc;%0.5; %new for full model
a2 = 360 + eS_c(8);%0.3;
alfa2 = 0 + eS_a(12);%0.12;
beta2 = 0 + eS_a(13);%0.16;
b2 = 0 + eS_c(9)*Kc;%0.14; %new for full model

d3 = 0 + eS_c(10);%0.2;
a3 = 0 + eS_c(11);%0.6;
alfa3 = 90 + eS_a(14);%0.11;
beta3 = 0 + eS_a(15)*Ka;%0.14; %new for full model
b3 = 0 + eS_c(12)*Kc;%0.2; %new for full model

d4 = 380 + eS_c(13);%0.2;
a4 = 0 + eS_c(14);%0.14;
alfa4 = -90 + eS_a(16);%0.13;
beta4 = 0 + eS_a(17)*Ka;%0.3; %new for full model
b4 = 0 + eS_c(15)*Kc;%0.13; %new for full model

d5 = 0 + eS_c(16);%0.15;
a5 = 0 + eS_c(17);%0.19;
alfa5 = 90 + eS_a(18);%0.14;
beta5 = 0 + eS_a(19)*Ka;%0.1; %new for full model
b5 = 0 + eS_c(18)*Kc;%0.11; %new for full model

d6 = 65 + eS_c(19);%0.4;
a6 = 0 + eS_c(20);%0.3;
alfa6 = 0 + eS_a(20);%0.12;
beta6 = 0 + eS_a(21)*Ka;%0.35; %new for full model
b6 = 0 + eS_c(21)*Kc;%0.18; %new for full model

M = [];

%% mean parameter mu and standard deviation parameter sigma.
rng(101994,'twister'); % random controller
ePmX = normrnd(0,0.15,[1,length(Th1)]); % random errors values instrument with normal distribution
rng(111960,'twister'); % random controller
ePmY = normrnd(0,0.15,[1,length(Th1)]); % random errors values instrument with normal distribution
rng(431951,'twister'); % random controller
ePmZ = normrnd(0,0.15,[1,length(Th1)]); % random errors values instrument with normal distribution

%% random errors values instrument with normal distribution
% rng(111960,'twister'); % random controller
% ePmX = 2*0.1*(randi([0 1],1,length(Th1)))-0.1; % mean parameter mu and standard deviation parameter sigma.
% rng(111960,'twister'); % random controller
% ePmY = 2*0.1*(randi([0 1],1,length(Th1)))-0.1; % mean parameter mu and standard deviation parameter sigma.
% rng(111960,'twister'); % random controller
% ePmZ = 2*0.1*(randi([0 1],1,length(Th1)))-0.1; % mean parameter mu and standard deviation parameter sigma.

for i=1:length(s1)
    
    RotZ0=[cosd(Th0) -sind(Th0) 0 0; sind(Th0) cosd(Th0) 0 0; 0 0 1 0; 0 0 0 1];
    Trans_d0=[1 0 0 0; 0 1 0 0; 0 0 1 d0; 0 0 0 1];
    Tran_a0=[1 0 0 a0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
    RotX0=[1 0 0 0; 0 cosd(alfa0) -sind(alfa0) 0; 0 sind(alfa0) cosd(alfa0) 0; 0 0 0 1];
    RotY0=[cosd(beta0) 0 sind(beta0) 0; 0 1 0 0; -sind(beta0) 0 cosd(beta0) 0;0 0 0 1];
    Trans_b0=[1 0 0 0; 0 1 0 b0; 0 0 1 0; 0 0 0 1];

    A0(:,:)=Tran_a0*Trans_b0*Trans_d0*RotZ0*RotY0*RotX0; %Full
%     A0(:,:,i)=RotZ0*Trans_d0*Tran_a0*RotX0*RotY0*Trans_b0;

    RotZ1(:,:)=[c1(i) -s1(i) 0 0; s1(i) c1(i) 0 0; 0 0 1 0; 0 0 0 1];
    Trans_d1=[1 0 0 0; 0 1 0 0; 0 0 1 d1; 0 0 0 1];
    Tran_a1=[1 0 0 a1; 0 1 0 0; 0 0 1 0; 0 0 0 1];
    RotX1=[1 0 0 0; 0 cosd(alfa1) -sind(alfa1) 0; 0 sind(alfa1) cosd(alfa1) 0; 0 0 0 1];
    RotY1=[cosd(beta1) 0 sind(beta1) 0; 0 1 0 0; -sind(beta1) 0 cosd(beta1) 0;0 0 0 1];
    Trans_b1=[1 0 0 0; 0 1 0 b1; 0 0 1 0; 0 0 0 1];
    
%     A1(:,:,i)=RotZ1(:,:,i)*Trans_d1*Tran_a1*RotX1;
    %A1=RotX1*Tran_a1*RotZ1*Trans_d1;
    A1(:,:)=RotZ1(:,:)*Trans_d1*Tran_a1*RotX1*RotY1*Trans_b1; %Full

    RotZ2(:,:)=[c2(i) -s2(i) 0 0; s2(i) c2(i) 0 0; 0 0 1 0; 0 0 0 1];
    Trans_d2=[1 0 0 0; 0 1 0 0; 0 0 1 d2; 0 0 0 1];
    Tran_a2=[1 0 0 a2; 0 1 0 0; 0 0 1 0; 0 0 0 1];
    RotX2=[1 0 0 0; 0 cosd(alfa2) -sind(alfa2) 0; 0 sind(alfa2) cosd(alfa2) 0; 0 0 0 1];
    RotY2=[cosd(beta2) 0 sind(beta2) 0; 0 1 0 0; -sind(beta2) 0 cosd(beta2) 0;0 0 0 1];
    Trans_b2=[1 0 0 0; 0 1 0 b2; 0 0 1 0; 0 0 0 1];

%     A2(:,:,i)=RotZ2(:,:,i)*Tran_a2*RotX2*RotY2; %Hayati
    % A2=RotZ2*Trans_d2*Tran_a2*RotX2; %DH
    %A2=RotX2*Tran_a2*RotZ2*Trans_d2;
    A2(:,:)=RotZ2(:,:)*Trans_d2*Tran_a2*RotX2*RotY2*Trans_b2; %Full

    RotZ3(:,:)=[c3(i) -s3(i) 0 0; s3(i) c3(i) 0 0; 0 0 1 0; 0 0 0 1];
    Trans_d3=[1 0 0 0; 0 1 0 0; 0 0 1 d3; 0 0 0 1];
    Tran_a3=[1 0 0 a3; 0 1 0 0; 0 0 1 0; 0 0 0 1];
    RotX3=[1 0 0 0; 0 cosd(alfa3) -sind(alfa3) 0; 0 sind(alfa3) cosd(alfa3) 0; 0 0 0 1];
    RotY3=[cosd(-90) 0 sind(-90) 0; 0 1 0 0; -sind(-90) 0 cosd(-90) 0;0 0 0 1];
    RotY3=[cosd(beta3) 0 sind(beta3) 0; 0 1 0 0; -sind(beta3) 0 cosd(beta3) 0;0 0 0 1];
    Trans_b3=[1 0 0 0; 0 1 0 b3; 0 0 1 0; 0 0 0 1];

    % A3=RotZ3*Tran_a3*RotX3*RotY3; %Hayati

%     A3(:,:,i)=RotZ3(:,:,i)*Trans_d3*Tran_a3*RotX3; %DH
    %A3=RotX3*Tran_a3*RotZ3*Trans_d3;
    A3(:,:)=RotZ3(:,:)*Trans_d3*Tran_a3*RotX3*RotY3*Trans_b3; %Full

    RotZ4(:,:)=[c4(i) -s4(i) 0 0; s4(i) c4(i) 0 0; 0 0 1 0; 0 0 0 1];
    Trans_d4=[1 0 0 0; 0 1 0 0; 0 0 1 d4; 0 0 0 1];
    Tran_a4=[1 0 0 a4; 0 1 0 0; 0 0 1 0; 0 0 0 1];
    RotX4=[1 0 0 0; 0 cosd(alfa4) -sind(alfa4) 0; 0 sind(alfa4) cosd(alfa4) 0; 0 0 0 1];
    RotY4=[cosd(beta4) 0 sind(beta4) 0; 0 1 0 0; -sind(beta4) 0 cosd(beta4) 0;0 0 0 1];
    Trans_b4=[1 0 0 0; 0 1 0 b4; 0 0 1 0; 0 0 0 1];

%     A4(:,:,i)=RotZ4(:,:,i)*Trans_d4*Tran_a4*RotX4;
    %A4=RotX4*Tran_a4*RotZ4*Trans_d4;
    A4(:,:)=RotZ4(:,:)*Trans_d4*Tran_a4*RotX4*RotY4*Trans_b4; %Full

    RotZ5(:,:)=[c5(i) -s5(i) 0 0; s5(i) c5(i) 0 0; 0 0 1 0; 0 0 0 1];
    Trans_d5=[1 0 0 0; 0 1 0 0; 0 0 1 d5; 0 0 0 1];
    Tran_a5=[1 0 0 a5; 0 1 0 0; 0 0 1 0; 0 0 0 1];
    RotX5=[1 0 0 0; 0 cosd(alfa5) -sind(alfa5) 0; 0 sind(alfa5) cosd(alfa5) 0; 0 0 0 1];
    RotY5=[cosd(beta5) 0 sind(beta5) 0; 0 1 0 0; -sind(beta5) 0 cosd(beta5) 0;0 0 0 1];
    Trans_b5=[1 0 0 0; 0 1 0 b5; 0 0 1 0; 0 0 0 1];

%     A5(:,:,i)=RotZ5(:,:,i)*Trans_d5*Tran_a5*RotX5;
    %A5=RotX5*Tran_a5*RotZ5*Trans_d5;
    A5(:,:)=RotZ5(:,:)*Trans_d5*Tran_a5*RotX5*RotY5*Trans_b5; %Full

    RotZ6(:,:)=[c6(i) -s6(i) 0 0; s6(i) c6(i) 0 0; 0 0 1 0; 0 0 0 1];
    Trans_d6=[1 0 0 0; 0 1 0 0; 0 0 1 d6; 0 0 0 1];
    Tran_a6=[1 0 0 a6; 0 1 0 0; 0 0 1 0; 0 0 0 1];
    RotX6=[1 0 0 0; 0 cosd(alfa6) -sind(alfa6) 0; 0 sind(alfa6) cosd(alfa6) 0; 0 0 0 1];
    RotY6=[cosd(beta6) 0 sind(beta6) 0; 0 1 0 0; -sind(beta6) 0 cosd(beta6) 0;0 0 0 1];
    Trans_b6=[1 0 0 0; 0 1 0 b6; 0 0 1 0; 0 0 0 1];

%     A6(:,:,i)=RotZ6(:,:,i)*Trans_d6*Tran_a6*RotX6;
    %A6=RotX6*Tran_a6*RotZ6*Trans_d6;
    A6(:,:)=RotZ6(:,:)*Trans_d6*Tran_a6*RotX6*RotY6*Trans_b6; %Full

    TI(:,:) = A0(:,:)*A1(:,:);
    TI(:,:) = TI(:,:)*A2(:,:);
    TI(:,:) = TI(:,:)*A3(:,:);
    TI(:,:) = TI(:,:)*A4(:,:);
    TI(:,:) = TI(:,:)*A5(:,:);
    TI(:,:) = TI(:,:)*A6(:,:);
    
    Mx(i) = round(TI(1,4),5) + ePmX(i);
    My(i) = round(TI(2,4),5) + ePmY(i);
    Mz(i) = round(TI(3,4),5) + ePmZ(i);
    
    
    M = [M Mx(i) My(i) Mz(i)];

end