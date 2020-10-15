clear all
close all
clc

%%
%Universidade de Brasilia
%LABROVIS
%Autor: Juan S. Toquica
%Ultima mod: outubro/2020
%Programa validar a calibração do IRB140 por sub-região
%%

tic

%% Leitura dataset
path = 'C:\Users\Juan\Documents\rovis_outubro_2020\kinematics_irb140';

file = 'IK-140_error_52_setembro.dat'; %controlar versão dos dataset
table_sr = readtable(fullfile(path, file));


it_s = 1;
%% Processo de validação da calibração

Th1=[];
% for it=[3 13 33 43]
for it=1:50
 %numero de SR - variavel dinamica
    
    row_end = it*35;
    row_ini = row_end - 34;
    
tab = table_sr{row_ini:row_end,1:9};
temp =1;
pos = 0;

for i=1:length(tab(:,1)) %i = 2, inicia em na posição 2 devido à posição (XYZ) da SR

    if(tab(i,4) > -180 && tab(i,4) < 180 && ...
            tab(i,5) > -90 && tab(i,5) < 110 && ...
            tab(i,6) > -230 && tab(i,6) < 50 &&...
            tab(i,7) > -200 && tab(i,7) < 200 &&...
            tab(i,8) > -115 && tab(i,8) < 115 &&...
            tab(i,9) > -400 && tab(i,9) < 400)
        pos(temp) = i;
        temp = temp + 1;
    end
end

% Em theta2 e theta3 foram incluidos os valores de -90 e 180
% para coincidir com os valores de XYZ do controlador IRC5 

% rng(18051989,'twister');
% cal_point = randi([2 length(pos)]);

x0 = tab(1,1);
y0 = tab(1,2);
z0 = tab(1,3);


xr=[];
yr=[];
zr=[];
r=[];

for i=2:length(pos)

    xr(i-1)=tab(i,1);
    yr(i-1)=tab(i,2);
    zr(i-1)=tab(i,3);
    
    r(i-1) = sqrt((xr(i-1)-x0)^2 + (yr(i-1)-y0)^2 + (zr(i-1)-z0)^2);

end

[val_r pos_r] = min(r);

Th1 = round(tab(pos,4),2); %-180:180
Th2 = round(tab(pos,5),2) - 90; %-90:110
Th3 = round(tab(pos,6),2) + 180; %-230:50
Th4 = round(tab(pos,7),2); %-200:200
Th5 = round(tab(pos,8),2); %-115:115
Th6 = round(tab(pos,9),2); %-400:400

% it_s = it_s + 1;

%Estanciar Theta
Theta = [Th1 Th2 Th3 Th4 Th5 Th6];



step = 1e-7;

[Mx My Mz M] = fk_irb140_med(Theta);
[Xn Yn Zn A Ap] = fk_irb140_base(Theta, step);
[X Y Z Ao Apo] = fk_irb140_base(Theta, 0);


%% Algoritmo LM
% Montar o Jacobiano

num_par = size(Ap);

J = [];
J = (Ao' - Ap') / step; 

[U,S,V] = svd(J);

deltaT = M - Ao;

error_ini = norm(deltaT);

I = eye(num_par(1),num_par(1));
k = 0.001;

% Primeira iteração do LM
pseudo_J = inv(J' * J + k*I) * J'; 

deltaX = zeros(6,1);

deltaX = deltaX - pseudo_J*deltaT';

[Xn Yn Zn A2 Ap2] = fk_irb140_base(Theta, step, deltaX);
[Xi Yi Zi Ao2 Apo2] = fk_irb140_base(Theta, 0, deltaX);

J2 = J;

deltaX2 = deltaX;

% Iterações do L-M

K0 = [];

for i=1:100

    deltaT2 = M - Ao2;

    error_new = norm(deltaT2);
    
    ctrl = sqrt((error_new-error_ini)'*(error_new-error_ini));
    
    if (error_new > error_ini)
        k = 1.05*k;
    else
        k = k/1.05;
    end

    pseudo_J2 = inv(J2' * J2 + k*I) * J2'; 

    deltaX2 = deltaX2 - pseudo_J2*deltaT2';
    
    [Xn Yn Zn A2 Ap2] = fk_irb140_base(Theta, step, deltaX2);
    [Xi Yi Zi Ao2 Apo2] = fk_irb140_base(Theta, 0, deltaX2);

    deltaTp = M - Ap2;
    
    J2 = (Ao2' - Ap2') / step; 
    
    [U2,S2,V2] = svd(J2);
    
    K2 = S2(1,1)/S2(num_par(1),num_par(1));

    sigma = S2;

    for n=1:num_par(1)
        sigma(n,n) = 0.001 / norm(J2(:,n)*J2(:,n)');
    end
    
    J0 = J2.*sigma;
    
    [U0,S0,V0] = svd(J0);
    
    K0(i) = S0(1,1)/S0(num_par(1),num_par(1));
    
    if(K0(i) < 100 && ctrl < error_new) %(Roth et al., 1987)
        [max_par(1), max_par(2)] = max(abs(V2(1:length(V2),length(V2))));
        max_par(3) = i;
        break;
    end
    
    error_ini = error_new;
    
end

%% Avaliação da calibração 
% Distancia Euclidiana

eM = [];
eMnew = [];
for i=1:length(Th1)

    eM(i) = sqrt((Mx(i)-X(i))^2 + (My(i)-Y(i))^2 + (Mz(i)-Z(i))^2);
    eMnew(i) = sqrt((Mx(i)-Xi(i))^2 + (My(i)-Yi(i))^2 + (Mz(i)-Zi(i))^2);
    eMnewX(i) = sqrt((Mx(i)-Xi(i))^2);
    eMnewY(i) = sqrt((My(i)-Yi(i))^2);
    eMnewZ(i) = sqrt((Mz(i)-Zi(i))^2);

end

desv = std(eM);
EmaxM = max(eM);
prom(it_s) = mean(eM);

desvNew = std(eMnew);
EmaxMnew = max(eMnew);
promNew(it_s) = mean(eMnew);

error_param(it_s,:) = deltaX2;
pos_sr(it_s,:) = [tab(1,1) tab(1,2) tab(1,3)];
% K_sr(it) = K0; 

data_plot = [desv prom(it_s) EmaxM; desvNew promNew(it_s) EmaxMnew]
bar(data_plot);

it_s = it_s + 1;

end


toc