function[X Y Z A Ap] = fk_irb140(Theta, step, delta_base)

delta = [-1.494164747
0.129063866
0.374533084
0.417463072
0.797890481
-0.344489435
-0.189015455
2.106442972
0.065325173
1.030348488
-0.391405204
-0.843142018
-0.097105952
-0.178059083
0.311048422
0.006167284
0.064785603
-0.425887308
-0.200792662
0.251968977
-0.120377199
];

    
if exist('delta_base')
    
    %Parâmetro Base-mundo
    Th0_i = 0 + delta_base(1); alfa0_i = 0 + delta_base(2);
    beta0_i = 0 + delta_base(3); a0_i = 0 + delta_base(4);
    b0_i = 0 + delta_base(5); d0_i = 0 + delta_base(6);
    
    %Parâmetros Junta_1
    Th1_i = Theta(:,1);% + delta(7);
    d1_i = 352;% + delta(8);
    a1_i = 70 + delta(1); alfa1_i = -90 + delta(2);

    %parâmetros Junta_2 - Hayati
    % d2 = 0;
    Th2_i = Theta(:,2) + delta(3); a2_i = 360 + delta(4);
    alfa2_i = 0 + delta(5); beta2_i = 0 + delta(6);

    %Parâmetros Junta_3
    Th3_i = Theta(:,3) + delta(7); d3_i = 0 + delta(8);
    a3_i = 0 + delta(9); alfa3_i = 90 + delta(10);

    %Parâmetros Junta_4
    Th4_i = Theta(:,4) + delta(11); d4_i = 380 + delta(12);
    a4_i = 0 + delta(13); alfa4_i = -90 + delta(14);

    %Parâmetros Junta_5
    Th5_i = Theta(:,5) + delta(15); d5_i = 0 + delta(16);
    a5_i = 0 + delta(17);
    alfa5_i = 90 + delta(18);

    %Parâmetros Junta_6
    Th6_i = Theta(:,6) + delta(19);
    d6_i = 65 + delta(20); 
    a6_i = 0 + delta(21); alfa6_i = 0;% + delta(28);
   
else
    
    %Parâmetro Base-mundo
    Th0_i = 0; alfa0_i = 0; beta0_i = 0;
    a0_i = 0; b0_i = 0; d0_i = 0;
   
    %Parâmetros Junta_1
    Th1_i = Theta(:,1);% + delta(7);
    d1_i = 352;% + delta(8);
    a1_i = 70 + delta(1); alfa1_i = -90 + delta(2);

    %parâmetros Junta_2 - Hayati
    % d2 = 0;
    Th2_i = Theta(:,2) + delta(3); a2_i = 360 + delta(4);
    alfa2_i = 0 + delta(5); beta2_i = 0 + delta(6);

    %Parâmetros Junta_3
    Th3_i = Theta(:,3) + delta(7); d3_i = 0 + delta(8);
    a3_i = 0 + delta(9); alfa3_i = 90 + delta(10);

    %Parâmetros Junta_4
    Th4_i = Theta(:,4) + delta(11); d4_i = 380 + delta(12);
    a4_i = 0 + delta(13); alfa4_i = -90 + delta(14);

    %Parâmetros Junta_5
    Th5_i = Theta(:,5) + delta(15); d5_i = 0 + delta(16);
    a5_i = 0 + delta(17);
    alfa5_i = 90 + delta(18);

    %Parâmetros Junta_6
    Th6_i = Theta(:,6) + delta(19);
    d6_i = 65 + delta(20); 
    a6_i = 0 + delta(21); alfa6_i = 0;% + delta(28);
    
end

if (step>0)
    
    param = zeros(6,1);
   
    for p=1:length(param)

        param(p) = param(p) + step;
        
        %Parâmetro Base-mundo
        Th0 = Th0_i + param(1); alfa0 = alfa0_i + param(2);
        beta0 = beta0_i + param(3); a0 = a0_i + param(4);
        b0 = b0_i + param(5); d0 = d0_i + param(6);
        
        %Parâmetros Junta_1
        Th1 = Th1_i;% + param(7);
        d1 = d1_i;% + param(8);
        a1 = a1_i; alfa1 = alfa1_i;

        %parâmetros Junta_2 - Hayati
        % d2 = 0;
        Th2 = Th2_i; a2 = a2_i;
        alfa2 = alfa2_i; beta2 = beta2_i;

        %Parâmetros Junta_3
        Th3 = Th3_i; d3 = d3_i;
        a3 = a3_i; alfa3 = alfa3_i;

        %Parâmetros Junta_4
        Th4 = Th4_i; d4 = d4_i;
        a4 = a4_i; alfa4 = alfa4_i;

        %Parâmetros Junta_5
        Th5 = Th5_i; d5 = d5_i;
        a5 = a5_i;
        alfa5 = alfa5_i;

        %Parâmetros Junta_6
        Th6 = Th6_i;
        d6 = d6_i;
        a6 = a6_i; alfa6 = alfa6_i;% + param(28);
        
        A = [];

        %% Calculo das Matrizes de Transformação Homogênea
        for i=1:length(Th1)
            
            %-----------------------------------------------------------------------------------
            RotZ0=[cosd(Th0) -sind(Th0) 0 0; sind(Th0) cosd(Th0) 0 0; 0 0 1 0; 0 0 0 1];
            Trans_d0=[1 0 0 0; 0 1 0 0; 0 0 1 d0; 0 0 0 1];
            Tran_a0=[1 0 0 a0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
            RotX0=[1 0 0 0; 0 cosd(alfa0) -sind(alfa0) 0; 0 sind(alfa0) cosd(alfa0) 0; 0 0 0 1];
            RotY0=[cosd(beta0) 0 sind(beta0) 0; 0 1 0 0; -sind(beta0) 0 cosd(beta0) 0;0 0 0 1];
            Trans_b0=[1 0 0 0; 0 1 0 b0; 0 0 1 0; 0 0 0 1];

            A0(:,:,i)=Tran_a0*Trans_b0*Trans_d0*RotZ0*RotY0*RotX0; %Full
%             A0(:,:,i)=RotZ0*Trans_d0*Tran_a0*RotX0*RotY0*Trans_b0;

            %-----------------------------------------------------------------------------------
            RotZ1(:,:,i)=[cosd(Th1(i)) -sind(Th1(i)) 0 0; sind(Th1(i)) cosd(Th1(i)) 0 0; 0 0 1 0; 0 0 0 1];
            Trans_d1=[1 0 0 0; 0 1 0 0; 0 0 1 d1; 0 0 0 1];
            Tran_a1=[1 0 0 a1; 0 1 0 0; 0 0 1 0; 0 0 0 1];
            RotX1=[1 0 0 0; 0 cosd(alfa1) -sind(alfa1) 0; 0 sind(alfa1) cosd(alfa1) 0; 0 0 0 1];

            A1(:,:,i)=RotZ1(:,:,i)*Trans_d1*Tran_a1*RotX1;
            %A1=RotX1*Tran_a1*RotZ1*Trans_d1;
            %-----------------------------------------------------------------------------------
            RotZ2(:,:,i)=[cosd(Th2(i)) -sind(Th2(i)) 0 0; sind(Th2(i)) cosd(Th2(i)) 0 0; 0 0 1 0; 0 0 0 1];
            %Trans_d2=[1 0 0 0; 0 1 0 0; 0 0 1 d2; 0 0 0 1];
            Tran_a2=[1 0 0 a2; 0 1 0 0; 0 0 1 0; 0 0 0 1];
            RotX2=[1 0 0 0; 0 cosd(alfa2) -sind(alfa2) 0; 0 sind(alfa2) cosd(alfa2) 0; 0 0 0 1];
            RotY2=[cosd(beta2) 0 sind(beta2) 0; 0 1 0 0; -sind(beta2) 0 cosd(beta2) 0;0 0 0 1];

            A2(:,:,i)=RotZ2(:,:,i)*Tran_a2*RotX2*RotY2; %Hayati
            % A2=RotZ2*Trans_d2*Tran_a2*RotX2; %DH
            %A2=RotX2*Tran_a2*RotZ2*Trans_d2;
            %-----------------------------------------------------------------------------------
            RotZ3(:,:,i)=[cosd(Th3(i)) -sind(Th3(i)) 0 0; sind(Th3(i)) cosd(Th3(i)) 0 0; 0 0 1 0; 0 0 0 1];
            Trans_d3=[1 0 0 0; 0 1 0 0; 0 0 1 d3; 0 0 0 1];
            Tran_a3=[1 0 0 a3; 0 1 0 0; 0 0 1 0; 0 0 0 1];
            RotX3=[1 0 0 0; 0 cosd(alfa3) -sind(alfa3) 0; 0 sind(alfa3) cosd(alfa3) 0; 0 0 0 1];
            RotY3=[cosd(-90) 0 sind(-90) 0; 0 1 0 0; -sind(-90) 0 cosd(-90) 0;0 0 0 1];

            % A3=RotZ3*Tran_a3*RotX3*RotY3; %Hayati

            A3(:,:,i)=RotZ3(:,:,i)*Trans_d3*Tran_a3*RotX3;
            %A3=RotX3*Tran_a3*RotZ3*Trans_d3;
            %-----------------------------------------------------------------------------------
            RotZ4(:,:,i)=[cosd(Th4(i)) -sind(Th4(i)) 0 0; sind(Th4(i)) cosd(Th4(i)) 0 0; 0 0 1 0; 0 0 0 1];
            Trans_d4=[1 0 0 0; 0 1 0 0; 0 0 1 d4; 0 0 0 1];
            Tran_a4=[1 0 0 a4; 0 1 0 0; 0 0 1 0; 0 0 0 1];
            RotX4=[1 0 0 0; 0 cosd(alfa4) -sind(alfa4) 0; 0 sind(alfa4) cosd(alfa4) 0; 0 0 0 1];

            A4(:,:,i)=RotZ4(:,:,i)*Trans_d4*Tran_a4*RotX4;
            %A4=RotX4*Tran_a4*RotZ4*Trans_d4;
            %-----------------------------------------------------------------------------------
            RotZ5(:,:,i)=[cosd(Th5(i)) -sind(Th5(i)) 0 0; sind(Th5(i)) cosd(Th5(i)) 0 0; 0 0 1 0; 0 0 0 1];
            Trans_d5=[1 0 0 0; 0 1 0 0; 0 0 1 d5; 0 0 0 1];
            Tran_a5=[1 0 0 a5; 0 1 0 0; 0 0 1 0; 0 0 0 1];
            RotX5=[1 0 0 0; 0 cosd(alfa5) -sind(alfa5) 0; 0 sind(alfa5) cosd(alfa5) 0; 0 0 0 1];

            A5(:,:,i)=RotZ5(:,:,i)*Trans_d5*Tran_a5*RotX5;
            %A5=RotX5*Tran_a5*RotZ5*Trans_d5;
            %-----------------------------------------------------------------------------------
            RotZ6(:,:,i)=[cosd(Th6(i)) -sind(Th6(i)) 0 0; sind(Th6(i)) cosd(Th6(i)) 0 0; 0 0 1 0; 0 0 0 1];
            Trans_d6=[1 0 0 0; 0 1 0 0; 0 0 1 d6; 0 0 0 1];
            Tran_a6=[1 0 0 a6; 0 1 0 0; 0 0 1 0; 0 0 0 1];
            RotX6=[1 0 0 0; 0 cosd(alfa6) -sind(alfa6) 0; 0 sind(alfa6) cosd(alfa6) 0; 0 0 0 1];

            A6(:,:,i)=RotZ6(:,:,i)*Trans_d6*Tran_a6*RotX6;
            %A6=RotX6*Tran_a6*RotZ6*Trans_d6;
            %-----------------------------------------------------------------------------------
            TI(:,:,i) = A0(:,:,i)*A1(:,:,i);
            TI(:,:,i) = TI(:,:,i)*A2(:,:,i);
            TI(:,:,i) = TI(:,:,i)*A3(:,:,i);
            TI(:,:,i) = TI(:,:,i)*A4(:,:,i);
            TI(:,:,i) = TI(:,:,i)*A5(:,:,i);
            TI(:,:,i) = TI(:,:,i)*A6(:,:,i);

            X(i) = TI(1,4,i);
            Y(i) = TI(2,4,i);
            Z(i) = TI(3,4,i);

            A = [A X(i) Y(i) Z(i)];

        end
        
        Ap(p,:) = A;

        param(p) = param(p) - step;
    
    end
    
else
    
    A = [];
    Ap = [];

    %% Calculo das Matrizes de Transformação Homogênea
    for i=1:length(Th1_i)


        %-----------------------------------------------------------------------------------
        RotZ0=[cosd(Th0_i) -sind(Th0_i) 0 0; sind(Th0_i) cosd(Th0_i) 0 0; 0 0 1 0; 0 0 0 1];
        Trans_d0=[1 0 0 0; 0 1 0 0; 0 0 1 d0_i; 0 0 0 1];
        Tran_a0=[1 0 0 a0_i; 0 1 0 0; 0 0 1 0; 0 0 0 1];
        RotX0=[1 0 0 0; 0 cosd(alfa0_i) -sind(alfa0_i) 0; 0 sind(alfa0_i) cosd(alfa0_i) 0; 0 0 0 1];
        RotY0=[cosd(beta0_i) 0 sind(beta0_i) 0; 0 1 0 0; -sind(beta0_i) 0 cosd(beta0_i) 0;0 0 0 1];
        Trans_b0=[1 0 0 0; 0 1 0 b0_i; 0 0 1 0; 0 0 0 1];

        A0(:,:,i)=Tran_a0*Trans_b0*Trans_d0*RotZ0*RotY0*RotX0; %Full
%         A0(:,:,i)=RotZ0*Trans_d0*Tran_a0*RotX0*RotY0*Trans_b0;
        
        %-----------------------------------------------------------------------------------
        RotZ1(:,:,i)=[cosd(Th1_i(i)) -sind(Th1_i(i)) 0 0; sind(Th1_i(i)) cosd(Th1_i(i)) 0 0; 0 0 1 0; 0 0 0 1];
        Trans_d1=[1 0 0 0; 0 1 0 0; 0 0 1 d1_i; 0 0 0 1];
        Tran_a1=[1 0 0 a1_i; 0 1 0 0; 0 0 1 0; 0 0 0 1];
        RotX1=[1 0 0 0; 0 cosd(alfa1_i) -sind(alfa1_i) 0; 0 sind(alfa1_i) cosd(alfa1_i) 0; 0 0 0 1];

        A1(:,:,i)=RotZ1(:,:,i)*Trans_d1*Tran_a1*RotX1;
        %A1=RotX1*Tran_a1*RotZ1*Trans_d1;
        %-----------------------------------------------------------------------------------
        RotZ2(:,:,i)=[cosd(Th2_i(i)) -sind(Th2_i(i)) 0 0; sind(Th2_i(i)) cosd(Th2_i(i)) 0 0; 0 0 1 0; 0 0 0 1];
        %Trans_d2=[1 0 0 0; 0 1 0 0; 0 0 1 d2; 0 0 0 1];
        Tran_a2=[1 0 0 a2_i; 0 1 0 0; 0 0 1 0; 0 0 0 1];
        RotX2=[1 0 0 0; 0 cosd(alfa2_i) -sind(alfa2_i) 0; 0 sind(alfa2_i) cosd(alfa2_i) 0; 0 0 0 1];
        RotY2=[cosd(beta2_i) 0 sind(beta2_i) 0; 0 1 0 0; -sind(beta2_i) 0 cosd(beta2_i) 0;0 0 0 1];

        A2(:,:,i)=RotZ2(:,:,i)*Tran_a2*RotX2*RotY2; %Hayati
        % A2=RotZ2*Trans_d2*Tran_a2*RotX2; %DH
        %A2=RotX2*Tran_a2*RotZ2*Trans_d2;
        %-----------------------------------------------------------------------------------
        RotZ3(:,:,i)=[cosd(Th3_i(i)) -sind(Th3_i(i)) 0 0; sind(Th3_i(i)) cosd(Th3_i(i)) 0 0; 0 0 1 0; 0 0 0 1];
        Trans_d3=[1 0 0 0; 0 1 0 0; 0 0 1 d3_i; 0 0 0 1];
        Tran_a3=[1 0 0 a3_i; 0 1 0 0; 0 0 1 0; 0 0 0 1];
        RotX3=[1 0 0 0; 0 cosd(alfa3_i) -sind(alfa3_i) 0; 0 sind(alfa3_i) cosd(alfa3_i) 0; 0 0 0 1];
        RotY3=[cosd(-90) 0 sind(-90) 0; 0 1 0 0; -sind(-90) 0 cosd(-90) 0;0 0 0 1];

        % A3=RotZ3*Tran_a3*RotX3*RotY3; %Hayati

        A3(:,:,i)=RotZ3(:,:,i)*Trans_d3*Tran_a3*RotX3;
        %A3=RotX3*Tran_a3*RotZ3*Trans_d3;
        %-----------------------------------------------------------------------------------
        RotZ4(:,:,i)=[cosd(Th4_i(i)) -sind(Th4_i(i)) 0 0; sind(Th4_i(i)) cosd(Th4_i(i)) 0 0; 0 0 1 0; 0 0 0 1];
        Trans_d4=[1 0 0 0; 0 1 0 0; 0 0 1 d4_i; 0 0 0 1];
        Tran_a4=[1 0 0 a4_i; 0 1 0 0; 0 0 1 0; 0 0 0 1];
        RotX4=[1 0 0 0; 0 cosd(alfa4_i) -sind(alfa4_i) 0; 0 sind(alfa4_i) cosd(alfa4_i) 0; 0 0 0 1];

        A4(:,:,i)=RotZ4(:,:,i)*Trans_d4*Tran_a4*RotX4;
        %A4=RotX4*Tran_a4*RotZ4*Trans_d4;
        %-----------------------------------------------------------------------------------
        RotZ5(:,:,i)=[cosd(Th5_i(i)) -sind(Th5_i(i)) 0 0; sind(Th5_i(i)) cosd(Th5_i(i)) 0 0; 0 0 1 0; 0 0 0 1];
        Trans_d5=[1 0 0 0; 0 1 0 0; 0 0 1 d5_i; 0 0 0 1];
        Tran_a5=[1 0 0 a5_i; 0 1 0 0; 0 0 1 0; 0 0 0 1];
        RotX5=[1 0 0 0; 0 cosd(alfa5_i) -sind(alfa5_i) 0; 0 sind(alfa5_i) cosd(alfa5_i) 0; 0 0 0 1];

        A5(:,:,i)=RotZ5(:,:,i)*Trans_d5*Tran_a5*RotX5;
        %A5=RotX5*Tran_a5*RotZ5*Trans_d5;
        %-----------------------------------------------------------------------------------
        RotZ6(:,:,i)=[cosd(Th6_i(i)) -sind(Th6_i(i)) 0 0; sind(Th6_i(i)) cosd(Th6_i(i)) 0 0; 0 0 1 0; 0 0 0 1];
        Trans_d6=[1 0 0 0; 0 1 0 0; 0 0 1 d6_i; 0 0 0 1];
        Tran_a6=[1 0 0 a6_i; 0 1 0 0; 0 0 1 0; 0 0 0 1];
        RotX6=[1 0 0 0; 0 cosd(alfa6_i) -sind(alfa6_i) 0; 0 sind(alfa6_i) cosd(alfa6_i) 0; 0 0 0 1];

        A6(:,:,i)=RotZ6(:,:,i)*Trans_d6*Tran_a6*RotX6;
        %A6=RotX6*Tran_a6*RotZ6*Trans_d6;
        %-----------------------------------------------------------------------------------
        TI(:,:,i) = A0(:,:,i)*A1(:,:,i);
        TI(:,:,i) = TI(:,:,i)*A2(:,:,i);
        TI(:,:,i) = TI(:,:,i)*A3(:,:,i);
        TI(:,:,i) = TI(:,:,i)*A4(:,:,i);
        TI(:,:,i) = TI(:,:,i)*A5(:,:,i);
        TI(:,:,i) = TI(:,:,i)*A6(:,:,i);

        X(i) = round(TI(1,4,i),5);
        Y(i) = round(TI(2,4,i),5);
        Z(i) = round(TI(3,4,i),5);

        A = [A X(i) Y(i) Z(i)];

    end
       
end

