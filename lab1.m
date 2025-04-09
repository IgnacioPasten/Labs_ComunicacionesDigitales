clc;
close all;
clear;


% Parámetros generales
fm = 100000;         % Frecuencia de muestreo muy alta (100 kHz)
tm = 1/fm;           % Periodo de muestreo
ls = 200;            % Número de muestras
t = (0:ls-1)*tm;     % Vector de tiempo para la señal original

% Señal senoidal m(t)
A = 1;               % Amplitud
f_c = 1000;          % Frecuencia de la señal senoidal (1 kHz)
m_t = A * sin(2*pi*f_c*t);  % Señal original

% Parámetros del muestreo PAM
f_s = 5000;          % Frecuencia de muestreo PAM (5 kHz)
t_s = 1/f_s;         % Periodo del muestreo
tau = 0.5 * t_s;     % Duración del pulso (50%)

% PAM Natural
pulse_train = zeros(1, length(t));
for n = 0:floor(t(end)/t_s)
    start_idx = find(t >= n*t_s, 1);
    end_idx = find(t >= n*t_s + tau, 1);
    if ~isempty(start_idx) && ~isempty(end_idx)
        pulse_train(start_idx:end_idx) = 1;
    end
end
pam_natural = m_t .* pulse_train;

% PAM Instantáneo (cuadrado)
pam_instant_rect = zeros(size(t));
for n = 0:floor(t(end)/t_s)
    sample_time = n * t_s;
    [~, idx_sample] = min(abs(t - sample_time));
    start_time = sample_time;
    end_time = sample_time + tau;
    idx_start = find(t >= start_time, 1);
    idx_end = find(t >= end_time, 1);
    if ~isempty(idx_start) && ~isempty(idx_end)
        pam_instant_rect(idx_start:idx_end) = m_t(idx_sample);
    end
end

%% Figura 1: Señales temporales (Subplots)
figure;

subplot(3,1,1);
plot(t*1000, m_t, 'b', 'LineWidth', 1.2);
title('Señal Original m(t)');
ylabel('Amplitud');
grid on;

subplot(3,1,2);
plot(t*1000, pam_natural, 'r', 'LineWidth', 1.2);
title('PAM Natural');
ylabel('Amplitud');
grid on;

subplot(3,1,3);
plot(t*1000, pam_instant_rect, 'g--', 'LineWidth', 1.2);
title('PAM Instantáneo (cuadrado)');
xlabel('Tiempo [ms]');
ylabel('Amplitud');
grid on;

sgtitle('Comparación Temporal de Señales');

%% Figura 2: Espectros (Subplots)
Nfft = 1024;
f = fm * (0:Nfft/2-1) / Nfft;

M_f = fft(m_t, Nfft);
PAM_nat_f = fft(pam_natural, Nfft);
PAM_inst_f = fft(pam_instant_rect, Nfft);

M_mag = abs(M_f(1:Nfft/2)) / max(abs(M_f));
PAM_nat_mag = abs(PAM_nat_f(1:Nfft/2)) / max(abs(PAM_nat_f));
PAM_inst_mag = abs(PAM_inst_f(1:Nfft/2)) / max(abs(PAM_inst_f));

figure;

subplot(3,1,1);
plot(f/1000, M_mag, 'b', 'LineWidth', 1.2);
title('Espectro de m(t)');
ylabel('Magnitud');
grid on;

subplot(3,1,2);
plot(f/1000, PAM_nat_mag, 'r', 'LineWidth', 1.2);
title('Espectro de PAM Natural');
ylabel('Magnitud');
grid on;

subplot(3,1,3);
plot(f/1000, PAM_inst_mag, 'g', 'LineWidth', 1.2);
title('Espectro de PAM Instantáneo');
xlabel('Frecuencia [kHz]');
ylabel('Magnitud');
grid on;

sgtitle('Espectros de Fourier');


N_bits = 6;                      
L = 2^N_bits;

m_t_norm = (m_t - min(m_t)) / (max(m_t) - min(m_t));
pam_inst_norm = (pam_instant_rect - min(pam_instant_rect)) / (max(pam_instant_rect) - min(pam_instant_rect));

min_val = min(pam_instant_rect);
max_val = max(pam_instant_rect);
q_step = (max_val - min_val) / L;
pam_pcm = floor((pam_instant_rect - min_val) / q_step) * q_step + min_val + q_step/2;
pam_pcm_norm = (pam_pcm - min(pam_pcm)) / (max(pam_pcm) - min(pam_pcm));

quant_error = pam_instant_rect - pam_pcm;

figure('Position', [100 100 1200 400]);

% =============================================
% Primer gráfico: Señal original, PAM y PCM
% =============================================
subplot(1,3,1);

plot(t*1000, m_t_norm, 'Color', [0 0.5 0.2], 'LineWidth', 2.2, 'DisplayName', 'Señal Original');
hold on;

plot(t*1000, pam_inst_norm, 'Color', [0.85 0.2 0.2], 'LineWidth', 1.6, 'DisplayName', 'PAM Instantáneo');

stairs(t*1000, pam_pcm_norm, 'Color', [0.1 0.3 0.9], 'LineWidth', 1.2, 'DisplayName', 'PCM (6 bits)');

title('Comparación de Señales', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Tiempo [ms]', 'FontSize', 10);
ylabel('Amplitud Normalizada', 'FontSize', 10);
legend('show', 'Location', 'northeast');
grid on;
xlim([0 max(t)*1000]);
ylim([-0.05 1.15]); % Margen ajustado
set(gca, 'FontSize', 10, 'GridAlpha', 0.3);


% =============================================
% Segundo gráfico: Señal PCM con stem - Versión Mejorada
% =============================================
subplot(1,3,2);

stem(t*1000, pam_pcm_norm, 'Color', [0 0.4 0.8], 'LineWidth', 1.2, ...
    'MarkerFaceColor', [0 0.4 0.8], 'MarkerEdgeColor', [0 0.2 0.6], ...
    'MarkerSize', 5, 'Marker', 'o');

title('Señal PCM Cuantizada (6 bits)', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Tiempo [ms]', 'FontSize', 10);
ylabel('Amplitud Normalizada', 'FontSize', 10);
grid on;
xlim([0 max(t)*1000]);
ylim([0 1.05]); 

hold on;
plot(xlim, [0 0], 'k:', 'LineWidth', 0.8);

sample_points = 1:round(t_s/tm):length(t);
plot(t(sample_points)*1000, pam_pcm_norm(sample_points), 'ro', ...
    'MarkerSize', 4, 'LineWidth', 1);

legend('Valores PCM', 'Puntos de muestreo', 'Location', 'southeast');

% =============================================
% Tercer gráfico: Error de cuantización
% =============================================
subplot(1,3,3);
plot(t*1000, quant_error, 'Color', [0.5 0.1 0.5], 'LineWidth', 1.8);
hold on;
fill([t*1000, fliplr(t*1000)], [quant_error, zeros(size(quant_error))], ...
     [0.8 0.6 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
title('Error de Cuantización', 'FontSize', 12);
xlabel('Tiempo [ms]', 'FontSize', 10);
ylabel('Error', 'FontSize', 10);
grid on;
xlim([0 max(t)*1000]);
ylim([-q_step q_step]*1.2);
set(gca, 'FontSize', 10);

plot(xlim, [0 0], 'k--', 'LineWidth', 1);


set(gcf, 'Position', [100 100 1200 350]);
set(gcf, 'Color', 'w');

figure('Position', [100 100 900 300]);
N_values = [4, 6, 8, 10];
colors = lines(length(N_values));

for i = 1:length(N_values)
    L_temp = 2^N_values(i);
    q_step_temp = (max_val - min_val)/L_temp;
    pam_pcm_temp = floor((pam_instant_rect - min_val)/q_step_temp) * q_step_temp + min_val + q_step_temp/2;
    quant_error_temp = pam_instant_rect - pam_pcm_temp;
    
    plot(t*1000, quant_error_temp, 'Color', colors(i,:), 'LineWidth', 1.5);
    hold on;
end

title('Error de Cuantización para distintos N', 'FontSize', 12);
xlabel('Tiempo [ms]');
ylabel('Error');
legend(cellstr(num2str(N_values')), 'Location', 'southeast');
grid on;
