%1)%

fc=1000;
fs= 100000;
tm=1/fs;
ls=200;
t=(0:ls-1)*tm;
m = sin(2* pi * fc * t);

%2)%

fs_2=5000;
Ts=1/fs_2;
ancho_pulso=0.5*Ts;
d=ancho_pulso/Ts;
pulsos_cuadrados=zeros(size(t));
n=0:Ts:max(t);
pam_natural=zeros(size(t));
for i = 1:length(n)
    pam_natural(t>=n(i) & t< n(i) + ancho_pulso)= m(t>=n(i) & t< n(i) + ancho_pulso);
end

%3)%

pam_instantaneo = zeros(size(t));
for i = 1:length(n)
    pam_instantaneo(t >= n(i) & t < n(i) + ancho_pulso) = sin(2* pi* fc* n(i));
end
figure;

% Primer gráfico: Señal senoidal
subplot(3,1,1);
plot(t, m);
title('Señal Senoidal');
xlabel('Tiempo (s)');
ylabel('Amplitud');
grid on;

% Segundo gráfico: Pulso cuadrático natural
subplot(3,1,2);
plot(t, pam_natural);
title('Pulso Cuadrático Natural');
xlabel('Tiempo (s)');
ylabel('Amplitud');
grid on;

% Tercer gráfico: Pulso instantáneo
subplot(3,1,3);
stem(t, pam_instantaneo);
title('Pulso Instantáneo');
xlabel('Tiempo (s)');
ylabel('Amplitud');
grid on;

sgtitle('Gráficos');
%Actividad%
%1)%
Nfft = length(t);                     
f = fs*(0:Nfft/2-1)/Nfft;            

M_f = fft(m);
PAM_nat_f = fft(pam_natural);
PAM_inst_f = fft(pam_instantaneo);

M_f = abs(M_f(1:Nfft/2));
PAM_nat_f = abs(PAM_nat_f(1:Nfft/2));
PAM_inst_f = abs(PAM_inst_f(1:Nfft/2));

% Primer gráfico: Señal senoidal
figure;
subplot(3,1,1);
plot(f, M_f/max(M_f));
title('Señal original');
xlabel('Tiempo (s)');
ylabel('Amplitud');
grid on;

% Segundo gráfico: Pulso cuadrático natural
subplot(3,1,2);
plot(f, PAM_nat_f/max(PAM_nat_f));
title('Forma natural');
xlabel('Tiempo (s)');
ylabel('Amplitud');
grid on;

% Tercer gráfico: Pulso instantáneo
subplot(3,1,3);
plot(f, PAM_inst_f/max(PAM_inst_f));
title('Forma instantanea');
xlabel('Tiempo (s)');
ylabel('Amplitud');
grid on;

sgtitle('Transformadas de Fourier');

%2)%
N=64;
L=2^N;
pcm_signal_inst = round((pam_instantaneo + 1) * (L - 1) / 2);
m_norm = (m - min(m)) / (max(m) - min(m));
pam_instantaneo_norm = (pam_instantaneo - min(pam_instantaneo)) / (max(pam_instantaneo) - min(pam_instantaneo));
pcm_signal_inst_norm = (pcm_signal_inst - min(pcm_signal_inst)) / (max(pcm_signal_inst) - min(pcm_signal_inst));
quantization_error_inst = pam_instantaneo_norm - ((2 * pcm_signal_inst / (L - 1)) - 1);

figure;
subplot(3,1,1);
hold on;
plot(t, m_norm);
plot(t, pam_instantaneo_norm);
plot(t, pcm_signal_inst_norm, '--');
title('PCM: Señal original, PAM instantáneo y PAM cuantificada');
legend('m(t)', 'PAM inst.', 'PCM');
xlabel('Tiempo (s)');
grid on;

% Gráfico del error
subplot(3,1,2);
plot(t, quantization_error_inst);
title('Error de cuantización (por muestra)');
xlabel('Tiempo (s)');
ylabel('Error');
grid on;
