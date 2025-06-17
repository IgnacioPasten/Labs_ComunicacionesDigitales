clear;
clc;
close all;

eb_n0_db = 1:11;

ber_bpsk = [0.0624, 0.0422, 0.0240, 0.0129, 0.0049, 0.0025, 0.0003, 0.0002, 0.00004, 0.000005, 0.0000004];
ber_qpsk = [0.0615, 0.0417, 0.0236, 0.0131, 0.0048, 0.0024, 0.0003, 0.0002, 0.00004, 0.000005, 0.0000004];
ber_8psk = [0.1862, 0.1479, 0.1096, 0.0759, 0.0468, 0.0257, 0.0132, 0.0065, 0.0023, 0.0008, 0.0002];

figure;

semilogy(eb_n0_db, ber_bpsk, '-o', 'LineWidth', 1.5, 'DisplayName', 'BPSK (Simulado)');
hold on; 
semilogy(eb_n0_db, ber_qpsk, '--s', 'LineWidth', 1.5, 'DisplayName', 'QPSK (Simulado)');
semilogy(eb_n0_db, ber_8psk, ':^', 'LineWidth', 1.5, 'DisplayName', '8-PSK (Simulado)');
hold off;


title('Rendimiento BER para Modulaciones PSK');
xlabel('Eb/N0 (dB)');
ylabel('Bit Error Rate (BER)');

legend('show');

grid on;

xlim([0, 12]);
ylim([1e-7, 1]);