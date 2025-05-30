clear all; close all; clc;

% Parámetros generales
Fs = 10000;
Rb = 100;
T_bit = 1/Rb;
N_bits = 10;
t1 = 0:1/Fs:N_bits*T_bit;

% Tren de pulsos
mt = square(2*pi*Rb*t1);
mt = (mt + 1)/2;
L = length(mt);
MT = fft(mt);
f_mt = (-L/2:L/2-1)*(Fs/L);
MT_mag = abs(fftshift(MT));

figure;
subplot(2,1,1);
plot(t1, mt); grid on;
title('Señal m(t) - Tren de pulsos');
xlabel('Tiempo (s)');
ylabel('Amplitud');

subplot(2,1,2);
plot(f_mt, MT_mag);
title('Transformada de Fourier |M(f)|');
xlabel('Frecuencia (Hz)');
ylabel('Magnitud');
grid on;

% Envolvente compleja g(t)
Ac = 1;
fs = 1000;
T = 1/fs;
t2 = 0:T:1-T;
m_t = sin(2*pi*5*t2) + 0.5*cos(2*pi*20*t2);
g_t = Ac * m_t;
G_f = fft(g_t);
f_g = linspace(-fs/2, fs/2, length(G_f));

figure;
plot(t2, m_t);
xlabel('Tiempo (s)');
ylabel('Amplitud');
title('Señal de Información m(t)');
grid on;

figure;
plot(f_g, abs(fftshift(G_f)));
xlabel('Frecuencia (Hz)');
ylabel('Magnitud');
title('Transformada de Fourier de g(t)');
grid on;

% Modulación FSK-4
Data = randi([0 1], 1, 10);
time = 0:9;
y = [];
j = 0;
for i = 1:2:length(Data)
    j = j + 1;
    if Data(i) == 0 && Data(i+1) == 0
        y(j) = 0;
    elseif Data(i) == 0 && Data(i+1) == 1
        y(j) = 1;
    elseif Data(i) == 1 && Data(i+1) == 1
        y(j) = 2;
    elseif Data(i) == 1 && Data(i+1) == 0
        y(j) = 3;
    end
end

figure;
stairs(time, Data, 'LineWidth', 2);
axis([0 9 -0.5 1.5]);
xlabel('Tiempo (bits)');
ylabel('Amplitud');
title(['Datos Binarios: ', num2str(Data)]);
grid on;

T_sym = 0:0.01:0.99;
f1 = sin(2 * pi * T_sym);
f2 = sin(8 * pi * T_sym);
f3 = sin(16 * pi * T_sym);
f4 = sin(32 * pi * T_sym);

signal = zeros(1, 100 * length(y));
for a = 1:length(y)
    idx = (a-1)*100 + 1 : a*100;
    switch y(a)
        case 0
            signal(idx) = f1;
        case 1
            signal(idx) = f2;
        case 2
            signal(idx) = f3;
        case 3
            signal(idx) = f4;
    end
end

t_fsk = 0:0.01:(length(signal)-1)*0.01;

figure;
plot(t_fsk, signal);
xlabel('Tiempo (s)');
ylabel('Amplitud');
title('Señal Modulada FSK-4');
grid on;

fs_FSK = 100;
G_FSK = fft(signal);
frequencies_FSK = linspace(-fs_FSK/2, fs_FSK/2, length(signal));

figure;
plot(frequencies_FSK, abs(fftshift(G_FSK)));
xlabel('Frecuencia (Hz)');
ylabel('Magnitud');
title('Transformada de Fourier de la Señal FSK');
grid on;
