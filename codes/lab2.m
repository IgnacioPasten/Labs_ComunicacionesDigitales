%CODE1 Respuesta al impulso y en frecuencia del filtro de coseno alzado

clear;
clc;
close all;

% Parámetros generales
f0 = 1; % Frecuencia central
B = 2 * f0; % Ancho de banda absoluto
t = linspace(0, 5, 1000); % Tiempo >= 0
f = linspace(-2*B, 2*B, 1000); % Frecuencia

% Valores de roll-off
alpha_values = [0, 0.25, 0.75, 1];

figure;
for i = 1:length(alpha_values)
alpha = alpha_values(i);
fDelta = alpha * f0;
f1 = f0 - fDelta;

% Respuesta en frecuencia He(f)
He = zeros(size(f));
for k = 1:length(f)
fk = abs(f(k));
if fk < f1
He(k) = 1;
elseif fk < f0 + fDelta
He(k) = 0.5 * (1 + cos(pi * (fk - f1) / (2 * fDelta)));
else
He(k) = 0;
end
end

% Respuesta al impulso he(t)
he = 2*f0 * (sinc(2*f0*t)) .* (cos(2*pi*fDelta*t) ./ (1 - (4*fDelta*t).^2));
he(isnan(he)) = 2*f0; 

% Gráficos
subplot(4, 2, 2*i-1);
plot(t, he);
title(['Respuesta al impulso - α = ', num2str(alpha)]);
xlabel('t'); ylabel('he(t)');

subplot(4, 2, 2*i);
plot(f, He);
title(['Respuesta en frecuencia - α = ', num2str(alpha)]);
xlabel('f'); ylabel('He(f)');
end


%CODE2 Solapamiento de pulsos para visualizar ISI

% Parámetros comunes
f0 = 1;                          
Ts = 1 / (2 * f0);              
span = 3;                       
t = linspace(-span*Ts*2, span*Ts*2, 1000); 
alpha_values = [0, 0.25, 0.75, 1]; 

% Figura general
figure;

for i = 1:length(alpha_values)
    alpha = alpha_values(i);

    % Función del pulso coseno alzado
    raised_cosine = @(t) 2*f0 .* sinc(2*f0*t) .* ...
        (cos(2*pi*alpha*f0*t) ./ (1 - (4*alpha*f0*t).^2));

    % Manejar posibles divisiones por cero
    t_eps = 1e-10;
    denom_zero = abs(1 - (4*alpha*f0*t).^2) < t_eps;
    rc = raised_cosine(t);
    rc(denom_zero) = 2*f0 * sinc(2*f0*t(denom_zero)); 

    % Suma de pulsos desplazados
    suma_pulsos = zeros(size(t));
    for k = -span:span
        suma_pulsos = suma_pulsos + raised_cosine(t - k*Ts);
    end

    % Subplot
    subplot(2, 2, i);
    plot(t, suma_pulsos, 'LineWidth', 1.5); hold on;
    yline(0, '--k');
    for k = -span:span
        xline(k*Ts, '--', 'Color', [0.6 0.6 0.6]);
    end
    title(['Solapamiento - \alpha = ', num2str(alpha)]);
    xlabel('Tiempo t'); ylabel('Amplitud');
    grid on;
end

sgtitle('Solapamiento de Pulsos Coseno Alzado para distintos valores de \alpha');

%CODE3 Visualizacion de Diagrama de Ojo

N = 10000;         % Número de bits a simular (10^4)
Ts = 1;            % Periodo de símbolo (s) - Normalizado
sps = 8;           % Muestras por símbolo (oversampling factor) 
alpha_values = [0, 0.25, 0.75, 1]; % Valores de roll-off
num_alpha = length(alpha_values);
SNR_dB = 20;       % Relación Señal a Ruido en dB 

% Filter span in symbols 
filter_span_symbols = 10;

% Generación de Datos
bits = randi([0 1], N, 1); % Generar N bits aleatorios

% Codificación NRZ-L (0 -> -1, 1 -> +1)
symbols_nrzl = 2*bits - 1;

% Simulación para cada alpha 
eyediagram_handles = gobjects(num_alpha, 1);

for i = 1:num_alpha
    alpha = alpha_values(i);

    current_alpha = alpha;
    if alpha == 0
        current_alpha = 1e-8;
        disp('Advertencia: Usando alpha cercano a 0 para rcosdesign.');
    end
    
    h_rrc = rcosdesign(current_alpha, filter_span_symbols, sps, 'normal');

    % Normalizar el filtro para tener ganancia unitaria en DC
    h_rrc = h_rrc / sum(h_rrc); 
    

    % Sobremuestrear los símbolos NRZ-L antes de filtrar
    symbols_upsampled = upsample(symbols_nrzl, sps);

    % Filtrar la señal sobremuestreada con el filtro coseno alzado
    tx_signal = filter(h_rrc, 1, symbols_upsampled);


    % Añadir Ruido Gaussiano Blanco Aditivo
    rx_signal = awgn(tx_signal, SNR_dB, 'measured'); 


    % Crear una nueva figura para cada diagrama de ojo
    eyediagram_handles(i) = figure;
    

    
    % Calcular el retardo introducido por el filtro FIR
    filter_delay = filter_span_symbols * sps / 2; 
    

    if length(rx_signal) > filter_delay
        signal_for_eye = rx_signal(filter_delay + 1 : end); 
    else
        warning('La señal filtrada es demasiado corta después de quitar el retardo del filtro.');
        signal_for_eye = rx_signal; 
    end

    % Graficar el diagrama de ojo mostrando 2 símbolos
    if ~isempty(signal_for_eye)
        eyediagram(signal_for_eye, 2*sps);
        
        title(sprintf('Diagrama de Ojo ($\\alpha = %.2f$, SNR = %d dB)', alpha, SNR_dB), 'Interpreter', 'latex');

        xlabel('Tiempo (dentro de 2T\_s)', 'Interpreter', 'latex'); 
        ylabel('Amplitud', 'Interpreter', 'latex');
        
        % Ajustar límites Y si es necesario para mejor visualización
        max_amp = max(abs(rx_signal));
         % Evitar error si max_amp es 0 o NaN
        if ~isnan(max_amp) && max_amp > 0 
            ylim([-1.5*max_amp, 1.5*max_amp]); 
        end
    else
        disp(['No se pudo generar el diagrama de ojo para alpha = ', num2str(alpha), ' debido a señal vacía.']);
    end

end

disp('Generación de diagramas de ojo completada.');


