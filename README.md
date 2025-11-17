# C-digo-tesis---Gurumendi-Carvajal-Christopher-Adonis
% Limpiar el espacio de trabajo y la ventana de comandos
clc;
clear all;
close all;
 
% --- Parámetros de la simulación ---
% Frecuencias de la portadora en Hz
fc_array = [2.5e9, 3.5e9, 4.5e9];
% Anchos de banda en Hz
B_array = [25e6, 100e6, 400e6];
% Velocidad del usuario en m/s (equivalente a 40 km/h)
velocidad_array = [11.11]; 
% Rango de distancias en metros
d = 10:10:1500;
% Velocidad de la luz en m/s
c = 3e8;
 
% --- DATOS AJUSTADOS PARA EL CANTÓN QUEVEDO ---
% Altura de la antena transmisora (Base Station) en metros
htx = 40; 
% Altura de la antena receptora (User Equipment) en metros
hrx = 1.5; 
% Altura media de los edificios
hr_buildings = 20; 
% Ancho medio de las calles
w_streets = 10; 
% Separación media entre edificios
b_buildings = 30; 
 
% Constante de Boltzmann y temperatura
k = 1.38e-23;
T = 290; % en Kelvin
% Factor de ruido del receptor
F = 10^(5/10); % 5 dB de factor de ruido (REALISTA)
% Potencia de transmisión en dBm
P_Tx_dBm = 40; % 40 dBm (10 W) (REALISTA)
 
% --- Parámetros de Simulación con Datos de Quevedo ---
% Desviación estándar para el efecto de shadowing
sigma_shadowing = 10;
% Distancia a la celda interferente (en metros)
d_interferer = 1500;
% Gap de la realidad para la capacidad (en dB)
gap_realidad_dB = 2;
 
% --- Eficiencia Espectral para Esquemas de Modulación Específicos (bits/s/Hz) ---
modulaciones = struct('QPSK', 2, 'QAM16', 4, 'QAM32', 5, 'QAM64', 6);
nombres_mod = {'QPSK', '16-QAM', '32-QAM', '64-QAM'};
 
% --- Bloque de almacenamiento de resultados para optimización ---
resultados_capacidad = [];
% Nombre actualizado del modelo COST-231
modelos_nombres = {'3GPP TR 38.901 (LoS/NLoS)', 'COST-231 Walfisch-ikegami', 'ITU-R P.1411', 'Xia-Bertoni'};
k_counter = 1;
 
% Colores para las líneas de modulación en todas las figuras
colores_mod_ac = {'g', 'b', 'c', 'm'}; 
 
% --- Bucle principal para cada combinación de parámetros ---
for i = 1:length(fc_array)
    fc = fc_array(i);
    
    for j = 1:length(B_array)
        B = B_array(j);
        
        % Potencia de ruido en dBm
        P_Ruido_W = k * T * B * F;
        P_Ruido_dBm = 10*log10(P_Ruido_W) + 30;
        
        % --- CÁLCULO DE LA PÉRDIDA DE TRAYECTORIA PARA CADA MODELO ---
        d_bp = 4 * (htx - 1) * (hrx - 1) * (fc/c);
        Lp_3gpp = zeros(size(d));
        for dist_idx = 1:length(d)
            if d(dist_idx) <= d_bp
                L_LoS = 22.0 * log10(d(dist_idx)) + 20 * log10(fc) - 28;
                Lp_3gpp(dist_idx) = L_LoS;
            else
                L_NLoS = 161.04 - 7.1*log10(b_buildings) + 7.5*log10(hr_buildings) - (24.37 - 3.7*(hr_buildings/htx)^2)*log10(htx) + 20*log10(fc/1e9) + 10*log10(d(dist_idx));
                Lp_3gpp(dist_idx) = L_NLoS;
            end
        end
        
        Lp_cost = 42.6 + 20 * log10(fc) + 26 * log10(d) - 147.56;
        k1 = 20;
        k2 = 10;
        k3 = 2;
        Lp_itu = k1*log10(d) + k2*log10(fc) + k3;
        h_eff = htx - hr_buildings;
        Lp_xia_bertoni = 42.6 + 20*log10(fc/1e9) + 26*log10(d) + 10*log10(h_eff);
        
        for l = 1:length(velocidad_array)
            velocidad_usuario = velocidad_array(l);
            factor_doppler_linear = 1 + (velocidad_usuario * fc / c)^2 / (B/2);
            modelos_celular = {Lp_3gpp, Lp_cost, Lp_itu, Lp_xia_bertoni};
            
            % --- Visualización en subgráficas de los 4 modelos ---
            % Se ajustó el tamaño de la figura aquí
            figure('Position', [100, 100, 1000, 600]); 
            sgtitle(['Capacidad del Canal (Frecuencia: ' num2str(fc/1e9) ' GHz, Ancho de Banda: ' num2str(B/1e6) ' MHz)']);
            
            for k_sub = 1:length(modelos_celular)
                Lp_actual = modelos_celular{k_sub};
                ganancia_fading = abs(randn(size(d)) + 1i*randn(size(d))) / sqrt(2);
                ganancia_fading_dB = 20*log10(ganancia_fading); 
                L_shadowing_dB = sigma_shadowing * randn(size(d));
                
                Lp_interferer_3gpp = 20 * log10(4*pi*d_interferer*fc/c) + 161.04 - 7.1*log10(b_buildings) + 7.5*log10(hr_buildings) - (24.37 - 3.7*(hr_buildings/htx)^2)*log10(htx) + 20*log10(fc/1e9); 
                P_interferer_dBm = P_Tx_dBm - Lp_interferer_3gpp;
                
                P_Rx_dBm = P_Tx_dBm - real(Lp_actual) - L_shadowing_dB + ganancia_fading_dB;
                P_Rx_W = 10.^(P_Rx_dBm/10) / 1000;
                P_interferer_W = 10.^(P_interferer_dBm/10) / 1000;
                
                SINR_ofdm_linear = P_Rx_W ./ (P_interferer_W + (10.^(P_Ruido_dBm/10)/1000) .* factor_doppler_linear);
                SINR_ofdm_dB = 10*log10(SINR_ofdm_linear);
                SINR_ofdm_dB_real = SINR_ofdm_dB - gap_realidad_dB;
                C_B_ofdm = log2(1 + 10.^(SINR_ofdm_dB_real./10));
                C_suavizada = smoothdata(C_B_ofdm, 'movmean', 50);
 
                % --- ALMACENAR RESULTADOS EN UNA ESTRUCTURA ---
                resultados_capacidad(k_counter).modelo = modelos_nombres{k_sub};
                resultados_capacidad(k_counter).fc = fc;
                resultados_capacidad(k_counter).B = B;
                resultados_capacidad(k_counter).capacidad = C_suavizada;
                resultados_capacidad(k_counter).promedio = mean(C_suavizada);
                k_counter = k_counter + 1;
                
                subplot(2, 2, k_sub);
                hold on;
                plot(d, C_B_ofdm, 'b', 'LineWidth', 1, 'DisplayName', 'Capacidad Real');
                plot(d, C_suavizada, 'k', 'LineWidth', 2.5, 'DisplayName', 'Tendencia de Capacidad');
                mod_fields = fieldnames(modulaciones);
                for m = 1:length(mod_fields)
                    mod_value = modulaciones.(mod_fields{m});
                    plot(d, ones(size(d)) * mod_value, 'Color', colores_mod_ac{m}, 'LineStyle', ':', 'LineWidth', 1.5, 'DisplayName', nombres_mod{m});
                end
                hold off;
                title(modelos_nombres{k_sub});
                xlabel('Distancia (m)');
                ylabel('Eficiencia Espectral (bits/s/Hz)');
                ylim([0 7]); 
                grid on;
                legend('show');
            end
        end
    end
end
 
% =========================================================================
% --- ANÁLISIS AUTOMÁTICO PARA ENCONTRAR LA MEJOR COMBINACIÓN ---
% =========================================================================
 
% Encontrar la combinación con la mayor capacidad promedio
[valor_max, indice_max] = max([resultados_capacidad.promedio]);
mejor_combinacion = resultados_capacidad(indice_max);
 
% Extraer los datos de la mejor combinación
mejor_modelo = mejor_combinacion.modelo;
mejor_fc = mejor_combinacion.fc;
mejor_B = mejor_combinacion.B;
mejor_capacidad_suavizada = mejor_combinacion.capacidad;
 
% --- Generar la figura extra con los cambios solicitados ---
% Se ajustó el tamaño de la figura aquí
figure('name', 'Canal Apropiado', 'Position', [100, 100, 1000, 600]); 
hold on;
plot(d, mejor_capacidad_suavizada, 'k', 'LineWidth', 2.5, 'DisplayName', 'Capacidad Real'); 
 
% Título actualizado aquí
title(sprintf('Canal Apropiado: %s (Frecuencia: %.1f GHz, Ancho de Banda: %.0f MHz)', mejor_modelo, mejor_fc/1e9, mejor_B/1e6)); 
xlabel('Distancia (m)');
ylabel('Eficiencia Espectral (bits/s/Hz)');
grid on;
 
% Graficar las líneas de modulación estáticas con sus colores
mod_fields = fieldnames(modulaciones);
for m = 1:length(mod_fields)
    mod_value = modulaciones.(mod_fields{m});
    plot(d, ones(size(d)) * mod_value, 'Color', colores_mod_ac{m}, 'LineStyle', '--', 'LineWidth', 1.0, 'DisplayName', nombres_mod{m});
end
 
legend('show', 'Location', 'best');
hold off;
