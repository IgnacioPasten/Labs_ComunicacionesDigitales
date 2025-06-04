clear all; close all; clc;

% --- SECCIÓN 1: Generación de Tren de Pulsos m(t) y su Transformada ---
disp('--- Sección 1: Generando Tren de Pulsos m(t) ---');
% Parámetros generales
Fs = 10000;     % Frecuencia de muestreo para la señal modulante
Rb = 100;       % Tasa de bits (bits por segundo)
T_bit = 1/Rb;   % Duración de un bit
N_bits = 10;    % Número de bits a simular
t1 = 0:1/Fs:(N_bits*T_bit - 1/Fs); % Vector de tiempo para la señal modulante

% Tren de pulsos (mensaje binario m(t))
mt = square(2*pi*Rb*t1/2, 50*(Fs/Rb)/Fs*100); % Se genera una onda cuadrada para simular los pulsos
mt(mt<0) = 0; 

L = length(mt);
MT = fft(mt);
f_mt = (-L/2:L/2-1)*(Fs/L);
MT_mag = abs(fftshift(MT));

figure;
subplot(2,1,1);
plot(t1, mt); grid on;
title('Señal m(t) - Tren de pulsos (0 o 1)');
xlabel('Tiempo (s)');
ylabel('Amplitud');
ylim([-0.2 1.2]);

subplot(2,1,2);
plot(f_mt, MT_mag);
title('Magnitud de la Transformada de Fourier |M(f)|');
xlabel('Frecuencia (Hz)');
ylabel('Magnitud');
grid on;

% --- SECCIÓN 2: Análisis de la Envolvente Compleja g(t) para OOK ---
disp('--- Sección 2: Análisis de la Envolvente Compleja g(t) para OOK ---');

% Parámetros para OOK
Ac_ook = 1; % Amplitud de la portadora para OOK

% Envolvente compleja g(t) para OOK
g_ook_t = Ac_ook * mt; 

% Transformada de Fourier de g_ook_t
L_g_ook = length(g_ook_t);
G_ook_f_mag = abs(fftshift(fft(g_ook_t)));
f_g_ook = (-L_g_ook/2:L_g_ook/2-1)*(Fs/L_g_ook); % Eje de frecuencias

% Gráficas para OOK g(t)
figure;
subplot(2,1,1);
plot(t1, g_ook_t);
title('Envolvente Compleja g_{OOK}(t) para OOK');
xlabel('Tiempo (s)');
ylabel('Amplitud');
grid on;
ylim([-0.1*Ac_ook Ac_ook+0.1*Ac_ook]); 

subplot(2,1,2);
plot(f_g_ook, G_ook_f_mag);
title('Magnitud de la Transformada de Fourier |G_{OOK}(f)| para OOK');
xlabel('Frecuencia (Hz)');
ylabel('Magnitud');
grid on;

% --- SECCIÓN 3: Análisis de la Envolvente Compleja g(t) para FSK BINARIA ---
disp('--- Sección 3: Análisis de la Envolvente Compleja g(t) para FSK Binaria ---');

% Parámetros para FSK Binaria
Ac_fsk = 1;     % Amplitud
delta_f = 200;  % Desviación de frecuencia en Hz (ej: 200 Hz).


% Señal modulante bipolar para FSK (a partir de mt de la Sección 1)
m_fsk_modulante = 2*mt - 1; 

% Cálculo de la fase phi(t) para g(t) = Ac * exp(j*phi(t))
phi_t_fsk = 2 * pi * delta_f * (cumsum(m_fsk_modulante) / Fs);

% Envolvente compleja g(t) para FSK Binaria
g_fsk_t = Ac_fsk * exp(1j * phi_t_fsk);

% Transformada de Fourier de g_fsk_t
L_g_fsk = length(g_fsk_t);
G_fsk_f_mag = abs(fftshift(fft(g_fsk_t)));
f_g_fsk = (-L_g_fsk/2:L_g_fsk/2-1)*(Fs/L_g_fsk); 

% Gráficas para FSK Binaria g(t)
figure;
subplot(3,1,1);
plot(t1, m_fsk_modulante);
title('Señal Modulante Bipolar m_{FSK}(t) para FSK (-1 o +1)');
xlabel('Tiempo (s)');
ylabel('Amplitud');
ylim([-1.2 1.2]);
grid on;

subplot(3,1,2);
plot(t1, angle(g_fsk_t)); % Graficamos la fase de g_fsk_t
title('Fase de la Envolvente Compleja \angle g_{FSK}(t) para FSK Binaria');
xlabel('Tiempo (s)');
ylabel('Fase (radianes)');
grid on;

subplot(3,1,3);
plot(f_g_fsk, G_fsk_f_mag);
title('Magnitud de la Transformada de Fourier |G_{FSK}(f)| para FSK Binaria');
xlabel('Frecuencia (Hz)');
ylabel('Magnitud');
grid on;

% --- SECCIÓN 4: Simulación de Modulación 4-FSK (del script original de tu compañero) ---
disp('--- Sección 4: Simulación de Modulación 4-FSK (Pasabanda) ---');

Data = randi([0 1], 1, 10); 
time_data_plot = 0:(length(Data)-1); 
y_symbols = []; 
j_symbol = 0;
for i_data = 1:2:length(Data)
    j_symbol = j_symbol + 1;
    if Data(i_data) == 0 && Data(i_data+1) == 0
        y_symbols(j_symbol) = 0;
    elseif Data(i_data) == 0 && Data(i_data+1) == 1
        y_symbols(j_symbol) = 1;
    elseif Data(i_data) == 1 && Data(i_data+1) == 1
        y_symbols(j_symbol) = 2;
    elseif Data(i_data) == 1 && Data(i_data+1) == 0
        y_symbols(j_symbol) = 3;
    end
end

figure;
stairs(time_data_plot, Data, 'LineWidth', 2);
axis([0 length(Data)-1 -0.5 1.5]);
xlabel('Índice de Bit');
ylabel('Valor del Bit');
title(['Datos Binarios Originales: ', num2str(Data)]);
grid on;

% Parámetros para la señal 4-FSK
samples_per_symbol = 100; % Muestras por cada símbolo FSK
T_symbol_duration_vector = (0:samples_per_symbol-1) / samples_per_symbol; 

% Frecuencias para cada símbolo (ejemplo, normalizadas respecto a una T_symbol_duration_vector de 0 a 1)
f_sym1 = sin(2 * pi * 2 * T_symbol_duration_vector); 
f_sym2 = sin(2 * pi * 4 * T_symbol_duration_vector);
f_sym3 = sin(2 * pi * 6 * T_symbol_duration_vector); 
f_sym4 = sin(2 * pi * 8 * T_symbol_duration_vector); 

signal_4fsk = zeros(1, samples_per_symbol * length(y_symbols));
for a_symbol = 1:length(y_symbols)
    idx_signal = (a_symbol-1)*samples_per_symbol + 1 : a_symbol*samples_per_symbol;
    switch y_symbols(a_symbol)
        case 0
            signal_4fsk(idx_signal) = f_sym1;
        case 1
            signal_4fsk(idx_signal) = f_sym2;
        case 2
            signal_4fsk(idx_signal) = f_sym3;
        case 3
            signal_4fsk(idx_signal) = f_sym4;
    end
end

% Eje de tiempo para la señal 4-FSK
dt_4fsk = 0.01 / samples_per_symbol;
t_4fsk = 0:dt_4fsk:(length(signal_4fsk)-1)*dt_4fsk;


figure;
plot(t_4fsk, signal_4fsk);
xlabel('Tiempo (s) (asumiendo T_{sym} = 0.01s)');
ylabel('Amplitud');
title('Señal Modulada 4-FSK (Pasabanda)');
grid on;

% Transformada de Fourier de la señal 4-FSK
Fs_4fsk_example = samples_per_symbol / 0.01; 
G_4FSK = fft(signal_4fsk);
frequencies_4FSK = linspace(-Fs_4fsk_example/2, Fs_4fsk_example/2, length(signal_4fsk));

figure;
plot(frequencies_4FSK, abs(fftshift(G_4FSK)));
xlabel('Frecuencia (Hz) (asumiendo Fs para 4-FSK = 10 kHz)');
ylabel('Magnitud');
title('Transformada de Fourier de la Señal 4-FSK (Pasabanda)');
grid on;

disp('--- Fin del Script ---');