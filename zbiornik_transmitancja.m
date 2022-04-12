% Model w postaci transmitancji (tf) ciągłej i dyskretnej,
% model w przestrzeni stanu (ss) ciągły i dyskretny,
% rysowanie odpowiedzi skokowych, badanie jakości dyskretyzacji

clear;

[punkt_pracy, parametry_modelu] = zbiornik_param;
[T_H,T_C,tau,tau_C,C0,alpha] = parametry_modelu();
[x10, x20, u10, u20, w10, w20]=punkt_pracy();

A = [1/(2*C0)*(-u10/(x10^2)-u20/(x10^2)-w10/(x10^2)+alpha/(2*sqrt(x10^3))), 0;
     1/C0*(-2*T_H*u10/(x10^3)-2*T_C*u20/(x10^3)-2*w20*w10/(x10^3)-alpha*(-3)*x20/(2*sqrt(x10^5))), -alpha/C0*(1/sqrt(x10^3))];
B = [1/(2*C0)*(1/x10), 1/(2*C0)*(1/x10);
     T_H/C0*(1/(x10^2)), T_C/C0*(1/(x10^2))];
C = [1, 0;
     0, 1];
D = [0, 0;
     0, 0];
 
const_h = 1/(2*C0)*((u10+u20+2*w10)/x10 - 3*alpha/(2*sqrt(x10)));
const_T = T_H/C0*2*u10/(x10^2) + T_C/C0*2*u20/(x10^2) + w10*w20*3/(C0*x10^2) - alpha*3*x20/(C0*2*sqrt(x10^3));

% Modele bez uwzględnienia opóźnień
SS0=ss(A,B,C,D);
TF0=tf(SS0);
[num, den] = tfdata(TF0);

% Uwzględnienie opóźnień w trasmitancji
TF = tf(num, den,'IODelay',[0, tau_C; tau, tau+tau_C]);

% Uwzględnienie opóźnień w modelu w przestrzeni stanu
DelayT(1) = struct('delay',120,'a',zeros(2),'b',[0, B(1,2); 0, B(2,2)],'c',zeros(2),'d',zeros(2));   % tauC=120
DelayT(2) = struct('delay',60,'a',zeros(2),'b',zeros(2),'c',[0,0;0,1],'d',zeros(2));  % tau=60
SS = delayss(A,[B(1,1), 0; B(2,1), 0],[1,0;0,0],zeros(2),DelayT);

Tp=1; % okres próbkowania
TFd = c2d(TF,Tp,'zoh'); % transmitancja dyskretna
opt = stepDataOptions('InputOffset',0,'StepAmplitude',1);
% figure; step(TF,opt); hold on; step(TFd,opt); grid on;
% legend('TF','TFd');

SSd = c2d(SS,Tp,'zoh'); % równania stanu dyskretne
figure; step(TFd,opt); hold on; step(SSd,opt); grid on;
legend('TFd','SSd');


% x0 = [54.39 34.05];
% t = 0:1:1000;
% u = zeros(length(t),2);
% u(1:200,1)=13;
% u(201:end,1)=13;
% u(1:200,2)=35;
% u(201:end,2)=35;
% lsim(TF,u,t,x0);
