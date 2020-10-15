clear all
close all
clc

%%
%Universidade de Brasilia
%LABROVIS
%Autor: Juan S. Toquica
%Ultima mod: outubro/2020
%Programa que faz a calibração do IRB140 segundo erros simulados com o
%conceito de sub-regiões. Por conjunto de sub-regiões.
%%

tic

%% Leitura dataset
path = 'C:\Users\Juan\Documents\rovis_outubro_2020\kinematics_irb140';

file = 'IK-140_error_52_setembro.dat'; %controlar versão dos dataset
table_sr = readtable(fullfile(path, file));


it_s = 1;
%% Processo de calibração

Th1 = [];
Th2 = [];
Th3 = [];
Th4 = [];
Th5 = [];
Th6 = [];

% for it=[2 7 42 47 72 77 92 97]
for it=[3 13 33 43]

    row_end = it*35;
    row_ini = row_end - 33;

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

it_s = it_s + 1;

end

%Estanciar Theta
Theta = [Th1 Th2 Th3 Th4 Th5 Th6];

step = 1e-7;

[Mx My Mz M] = fk_irb140_med(Theta);
[Xn Yn Zn A Ap] = fk_irb140(Theta, step);
[X Y Z Ao Apo] = fk_irb140(Theta, 0);


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

% Primeira iteração do L-M
pseudo_J = inv(J' * J + k*I) * J'; 

deltaX = zeros(num_par(1),1);

deltaX = deltaX - pseudo_J*deltaT';

[Xn Yn Zn A2 Ap2] = fk_irb140(Theta, step, deltaX);
[Xi Yi Zi Ao2 Apo2] = fk_irb140(Theta, 0, deltaX);

J2 = J;

deltaX2 = deltaX;

%Iterações do L-M

for i=1:100

    deltaT2 = M - Ao2;

    error_new = norm(deltaT2);
    
    ctrl = sqrt((error_new-error_ini)'*(error_new-error_ini));
    
    if (error_new > error_ini)
        k = 10*k;
    else
        k = k/10;
    end

    pseudo_J2 = inv(J2' * J2 + k*I) * J2'; 

    deltaX2 = deltaX2 - pseudo_J2*deltaT2';
    
    [Xn Yn Zn A2 Ap2] = fk_irb140(Theta, step, deltaX2);
    [Xi Yi Zi Ao2 Apo2] = fk_irb140(Theta, 0, deltaX2);

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
    
    if(K0(i) < 100 && ctrl < error_new)
        [max_par(1), max_par(2)] = max(abs(V2(7:length(V2),length(V2))));
        max_par(3) = i;
        break;
    end
    
    error_ini = error_new;
    
end

eM = 0;
for i=1:length(Th1)

    eM(i) = sqrt((Mx(i)-X(i))^2 + (My(i)-Y(i))^2 + (Mz(i)-Z(i))^2);
    eMnew(i) = sqrt((Mx(i)-Xi(i))^2 + (My(i)-Yi(i))^2 + (Mz(i)-Zi(i))^2);
    
    eMx(i) = sqrt((Mx(i)-X(i))^2);
    eMy(i) = sqrt((My(i)-Y(i))^2);
    eMz(i) = sqrt((Mz(i)-Z(i))^2);
    
    eMnewX(i) = sqrt((Mx(i)-Xi(i))^2);
    eMnewY(i) = sqrt((My(i)-Yi(i))^2);
    eMnewZ(i) = sqrt((Mz(i)-Zi(i))^2);

end

desv = std(eM);
EmaxM = max(eM);
prom = mean(eM);

desvNew = std(eMnew);
EmaxMnew = max(eMnew);
promNew = mean(eMnew);
promNewX = mean(eMnewX);
promNewY = mean(eMnewY);
promNewZ = mean(eMnewZ);

promX = mean(eMx);
promY = mean(eMy);
promZ = mean(eMz);

error_param = deltaX2;
pos_sr = [tab(1,1) tab(1,2) tab(1,3)];

data_plot = [desv prom EmaxM; desvNew promNew EmaxMnew]
bar(data_plot);

%% Gerar tabela com os parametros de erro 27/abril/2020

%%Perguntar se salvar dados ou nao


prompt = 'Salvar dataset? (y/n)';
answer_u = input(prompt,'s');


if (answer_u == 'y')
    
    a=datetime('now','Format','yyyy-MM-dd''T''HHmmss');
    a_data = datestr(a,'yyyymmddTHHMMSS');
    file = sprintf('SR_error_parameters_%s.dat',a_data);
    
    %XYZ de cada subregião
    
    pos_sr = table2array(pos_sr);
    
    if (num_par(1) == 27) %só 23 parâmetros melhorando o condicionamento
        
        TableSR1 = table(pos_sr(:,1), pos_sr(:,2), pos_sr(:,3), error_param(1),error_param(2),error_param(3),error_param(4),...
                    error_param(5),error_param(6),error_param(7),error_param(8),...
                    error_param(9),error_param(10),error_param(11),error_param(12),...
                    error_param(13),error_param(14),error_param(15),error_param(16),...
                    error_param(17),error_param(18),error_param(19),error_param(20),...
                    error_param(21),error_param(22),error_param(23),error_param(24),...
                    error_param(25),error_param(26),error_param(27),...
        'VariableNames',{'X' 'Y' 'Z' 'delta_Th0' 'delta_alfa0' 'delta_beta0' 'delta_a0' 'delta_b0' 'delta_d0'...
                        'delta_a1' 'delta_alfa1' 'delta_Th2' 'delta_a2' 'delta_alfa2' 'delta_beta2'...
                        'delta_Th3' 'delta_d3' 'delta_a3' 'delta_alfa3' 'delta_Th4' 'delta_d4'...
                        'delta_a4' 'delta_alfa4' 'delta_Th5' 'delta_d5' 'delta_a5' 'delta_alfa5' 'delta_Th6'...
                        'delta_d6' 'delta_a6'});  
        
    else %todos os parametros
       
        TableSR1 = table(pos_sr(:,1), pos_sr(:,2), pos_sr(:,3), error_param(1),error_param(2),error_param(3),error_param(4),...
                    error_param(5),error_param(6),error_param(7),error_param(8),...
                    error_param(9),error_param(10),error_param(11),error_param(12),...
                    error_param(13),error_param(14),error_param(15),error_param(16),...
                    error_param(17),error_param(18),error_param(19),error_param(20),...
                    error_param(21),error_param(22),error_param(23),error_param(24),...
        'VariableNames',{'X' 'Y' 'Z' 'delta_Th1' 'delta_d1' 'delta_a1' 'delta_alfa1' 'delta_Th2' 'delta_d2' 'delta_alfa2'...
                        'delta_beta2' 'delta_Th3' 'delta_d3' 'delta_a3' 'delta_alfa3' 'delta_Th4' 'delta_d4'...
                        'delta_a4' 'delta_alfa4' 'delta_Th5' 'delta_d5' 'delta_a5' 'delta_alfa5' 'delta_Th6'...
                        'delta_d6' 'delta_a6' 'delta_alfa6'}); 
    end
                  
    writetable(TableSR1,file); %


else

    answer_u = 'n';
    
end

toc