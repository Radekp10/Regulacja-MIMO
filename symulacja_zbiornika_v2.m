% Skrypt do symulacji zbiornika - badanie wyjścia dla zadanej trajektorii
% sterowań i zakłóceń

clear;

% Pobranie uchwytu do funkcji symulującej zbiornik
[punkt_pracy, parametry_modelu]=zbiornik_param;
[symulacja_obiektu]=zbiornik_sym;

% Nominalny punkt pracy
[x10, x20, u10, u20, w10, w20, y10, y20]=punkt_pracy();

% Parametry modelu
[T_H,T_C,tau,tau_C,C,alpha] = parametry_modelu();

n=2e3; % liczba kroków h
h=1; % wielkość kroku [s]

czy_zlin=true; % czy model liniowy czy zlinearyzowany

% Osie do wykresu
tu=0:1:(n-1);
ty=1:1:n;

x1=zeros(n,1);
x2=zeros(n,1);
y1=zeros(n,1);
y2=zeros(n,1);
u1=zeros(n-1,1);
u2=zeros(n-1,1);
w1=zeros(n-1,1);
w2=zeros(n-1,1);

% Początkowy punkt pracy
x1(1)=x10;
x2(1)=x20;
y1(1)=y10;
y2(1)=y20;

% Zadana trajektoria sterowań i zakłóceń
u1(1:n-1)=u10;
% u1(1000:end)=25;
u2(1:n-1)=u20;
w1(1:n-1)=w10;
w2(1:n-1)=w20;


% Główna pętla symulacji obiektu    
for i=2:n
    if i-1-tau_C > 0
        u2_opoz=u2(i-1-tau_C);
    else
        u2_opoz=u2(1);
    end
    
    if i-1-tau > 0
        x2_opoz=x2(i-1-tau);
    else
        x2_opoz=x2(1);
    end

    [y, x_nast] = symulacja_obiektu(u1(i-1), u2_opoz, w1(i-1), w2(i-1), x1(i-1), x2(i-1), x2_opoz, h, czy_zlin);
    x1(i)=x_nast(1);
    x2(i)=x_nast(2);
    y1(i)=y(1);
    y2(i)=y(2);

end    
 

% Wykresy
f=figure;
if czy_zlin 
    f.Name='Symulacja z modelem zlinearyzowanym';
else
    f.Name='Symulacja z modelem nieliniowym';
end
subplot(3,2,1)
stairs(tu, [u10; u1])
grid on;
title("Sygnal wejsciowy $u_1$",'interpreter','latex');
xlabel("Numer kroku",'interpreter','latex');
ylabel("$u_1$ $[cm^{3}/s]$",'interpreter','latex');

subplot(3,2,2)
stairs(tu, [u20; u2])
grid on;
title("Sygnal wejsciowy $u_2$",'interpreter','latex');
xlabel("Numer kroku",'interpreter','latex');
ylabel("$u_2$ $[cm^{3}/s]$",'interpreter','latex');

subplot(3,2,3);
plot(ty, [y1])
grid on;
title("Sygnal wyjsciowy $y_1$",'interpreter','latex');
xlabel("Numer kroku",'interpreter','latex');
ylabel("$y_1$ $[cm]$",'interpreter','latex');
% ylim([0.9*y1(1); 1.1*y1(1)]);

subplot(3,2,4);
plot(ty, [y2])
grid on;
title("Sygnal wyjsciowy $y_2$",'interpreter','latex');
xlabel("Numer kroku",'interpreter','latex');
ylabel("$y_2$ $[^\circ C]$",'interpreter','latex');
% ylim([0.9*y2(1); 1.1*y2(1)]);

subplot(3,2,5);
plot(ty, [x1])
grid on;
title("Zmienna stanu $x_1$",'interpreter','latex');
xlabel("Numer kroku",'interpreter','latex');
ylabel("$x_1$ $[cm]$",'interpreter','latex');
% ylim([0.9*x1(1); 1.1*x1(1)]);

subplot(3,2,6);
plot(ty, [x2])
grid on;
title("Zmienna stanu $x_2$",'interpreter','latex');
xlabel("Numer kroku",'interpreter','latex');
ylabel("$x_2$ $[^\circ C]$",'interpreter','latex');
% ylim([0.9*x2(1); 1.1*x2(1)]);
