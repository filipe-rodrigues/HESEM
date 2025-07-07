%% HESEM: HOLOMORPHIC EMBEDDING STATE ESTIMATION METHOD USING PHASOR MEASUREMENT UNITS

clear all ; close all
warning off
format short
clc

%% Initialization

tol = 1e-3;   % Tolerance
rec_max = 50; % Maximum number of Recursions
MC = 100;      % For Monte Carlo simulation
on = 1;       % Introduce aleatory measurement errors (yes(1)/no(2))
tic
%% Case studies

Sistema_estudado = 14; % 14 for the 14-bus; 33 for the 33-bus; 107 for the 107-bus
Cenario_estudado = 2;  % Scenario 1 or 2 accordng to the article

if Sistema_estudado == 14
    Dados_14_bus;
elseif Sistema_estudado == 33
    Dados_33_bus;
elseif Sistema_estudado == 107
    Dados_107_bus;
end

%% Power flow

[V,Skm,Sk,Ikm,Ybarra,NB,bbus,ykm,TAP] = HELM(tol,rec_max,DLIN,DBAR,Sb,Zbase); % power flow solution
V_cal = V; % actual values

%% HESEM

tic

Tensao_estimada = zeros(NB,MC);
contador = zeros(1,MC);

for iter = 1:MC
    
    %% Measurements

    medidas = Medicoes(Sistema_estudado,Cenario_estudado);

    Q_medidas = length(medidas(:,1)); % number of measurements

    tipo = medidas(:,2);              % type of measurement

    de = medidas(:,3);                % from bus
    para = medidas(:,4);              % to bus
    Rii = eye(length(medidas(:,5)));

    z = zeros(Q_medidas,1);

    for i = 1:Q_medidas

        if tipo(i) == 1     
            z(i) = V(de(i)) * (1+0.0004*randn*on);           % Vk measurement

        elseif tipo(i) == 2 
            z(i) = Sk(de(i)) * (1+0.001*randn*on);           % Sk measurement

        elseif tipo(i) == 3 
            z(i) = Skm(de(i),para(i)) * (1+0.001*randn*on);  % Skm measurement

        else                
            z(i) = Ikm(de(i),para(i)) * (1+0.0008*randn*on); % Ikm measurement
        end

    end

    %% State estimation - HESEM

    % Matrix H:

    Ysh = sum(Ybarra.').';  % Ybus -only shunt elements

    Ytr = Ybarra-diag(Ysh); % Ybus without shunt elements

    Y_rec = zeros(Q_medidas,NB);

    for i = 1:Q_medidas

        if tipo(i) == 1
            Y_rec(i, de(i))   = 1;

        elseif tipo(i) == 2
            Y_rec(i,:) = Ytr(de(i),:);

        elseif tipo(i) == 3
            Y_rec(i, de(i))     =   ykm(de(i),para(i)) * (1/TAP(de(i),para(i)))^2 + 1j*bbus(de(i),para(i));
            Y_rec(i, para(i))   =  -ykm(de(i),para(i)) * (1/TAP(de(i),para(i)));

        else
            Y_rec(i, de(i))     =   ykm(de(i),para(i)) * (1/TAP(de(i),para(i)))^2 + 1j*bbus(de(i),para(i));
            Y_rec(i, para(i))   =  -ykm(de(i),para(i)) * (1/TAP(de(i),para(i)));

        end
    end

    Wii=inv(Rii);

    Y_rec = inv((transpose(Y_rec)*Wii*Y_rec))*transpose(Y_rec)*Wii;

    % Wk calculation:

    coeficientes_W = @(v,w) -diag(v(:,2:end)*fliplr(w).'); 

    Vh_estimado           = ones(NB,1);
    W_estimado            = ones(NB,1);
    k            = 2;
    delS         = 1;

    V_old_estimado = ones(NB,1);

    V_est = zeros(NB,1);

    while (delS > tol)&&(k < rec_max)

        delta_z = calcula_Delta_Z(k,z,W_estimado,medidas,de,Ysh,Vh_estimado,ykm,TAP,bbus,para);

        Vh_estimado(:,k)  = Y_rec*delta_z;
        W_estimado(:,k)   = coeficientes_W(Vh_estimado,W_estimado);

        k        = k+1;

        for t = 1:length(Vh_estimado(:,1))
            V_est(t,1) = Aproximantes_Pade(Vh_estimado(t,:));
        end

        delS = max(abs(V_est-V_old_estimado));

        V_old_estimado = V_est;
    end

    Tensao_estimada(:,iter) = V_est;
    contador(iter) = k;
end

time = toc;

% For the figures
set(0,'defaultAxesFontName', 'times');
set(0,'defaultTextFontName', 'times');
set(0,'defaultAxesFontSize',14);
set(0,'defaultTextFontSize',14);

fprintf('Computer simulation time: %.2f\n', time)
fprintf('Number of recursions - HESEM: %.2f\n', mean(contador))

V_calculado = V_cal;
V_estimado = mean(Tensao_estimada, 2);  % vetor estimado médio

% Magnitude percentage error
Erro_V = abs(abs(V_calculado) - abs(V_estimado)) ./ abs(V_calculado) * 100;

% Angle asolute error
Erro_Angulo = abs(rad2deg(angle(V_calculado)) - rad2deg(angle(V_estimado)));

% Exibir resultados
fprintf('Magnitude percentage error: %.3f%%\n', max(Erro_V))
fprintf('Angle absolute error: %.3f graus\n', max(Erro_Angulo))


%--------------------------------------------------------------------------
% Plots
%--------------------------------------------------------------------------

my_blue = [0 0 .5];
my_red = [1 0 0];
barras = linspace(1,NB,NB);

% Figure 1: Voltage Magnitude Comparison
figure(1)
if Sistema_estudado == 107
    set(gcf, 'Position', [100, 100, 700, 400]); % [x_pos, y_pos, width, height]
end

plot(barras, abs(V_calculado), 'o-', 'LineWidth', 1.2, 'Color', my_red, 'MarkerFaceColor', my_red, 'MarkerSize', 4.5); 
hold on
plot(barras, abs(V_estimado), '--x', 'LineWidth', 1.2, 'Color', my_blue, 'MarkerSize', 4.5);
grid on
legend('Calculated Voltage', 'Estimated Voltage', 'Location', 'best') 
xlabel('Buses')
ylabel('V (pu)')
xlim([1,NB])
% title('Comparison between Calculated and Estimated Voltage')
hold off

% Figure 2: Voltage Angle Comparison
figure(2)
if Sistema_estudado == 107
    set(gcf, 'Position', [100, 100, 700, 400]); % [x_pos, y_pos, width, height]
end

plot(barras, rad2deg(angle(V_calculado)), 'o-', 'LineWidth', 1.2, 'Color', my_red, 'MarkerFaceColor', my_red, 'MarkerSize', 4.5); 
hold on
plot(barras, rad2deg(angle(V_estimado)), '--x', 'LineWidth', 1.2, 'Color', my_blue, 'MarkerSize', 4.5); 
grid on
legend('Calculated Angle', 'Estimated Angle', 'Location', 'best')
xlabel('Buses')
ylabel('Angle (º)')
xlim([1,NB])
% title('Comparison between Calculated and Estimated Angles')
hold off


% Figure 3: Voltage Magnitude Absolute Error
figure(3)
if Sistema_estudado == 107
    set(gcf, 'Position', [100, 100, 700, 400]); % [x_pos, y_pos, width, height]
end

plot(barras, abs(Erro_V), 'LineWidth', 1.2, 'Color', 'black', ...
    'Marker', 'd', 'MarkerFaceColor', 'black', 'MarkerSize', 4);
grid on
xlabel('Buses')
ylabel('Absolute Voltage Error (%)')
% title('Absolute Voltage Error (PMUs otimamente alocadas)')
xlim([1, NB])
hold off

% Figure 4: Voltage Angle Absolute Error
figure(4)
if Sistema_estudado == 107
    set(gcf, 'Position', [100, 100, 700, 400]); % [x_pos, y_pos, width, height]
end

plot(barras, abs(Erro_Angulo), 'LineWidth', 1.2, 'Color', 'black', ...
    'Marker', 'd', 'MarkerFaceColor', 'black', 'MarkerSize', 4);
grid on
xlabel('Buses')
ylabel('Absolute Angle Error (°)')
% title('Absolute Angle Error (PMUs otimamente alocadas)')
xlim([1, NB])
hold off
