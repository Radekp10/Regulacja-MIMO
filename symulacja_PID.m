% Symulacja obiektu z regulatorem dwupętlowym PID

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

czy_zlin=false; % czy model liniowy czy zlinearyzowany
czy_odsprz=false; % czy PID z odsprzęganiem czy bez

% Osie do wykresu
tu=0:1:(n-1);
ty=1:1:n;

x1=zeros(n,1);
x2=zeros(n,1);
y1=zeros(n,1);
y2=zeros(n,1);
y1zad=zeros(n,1);
y2zad=zeros(n,1);
u1=zeros(n-1,1);
u2=zeros(n-1,1);
uR1=zeros(n-1,1);
uR2=zeros(n-1,1);
w1=zeros(n-1,1);
w2=zeros(n-1,1);

% Początkowy punkt pracy
x1(1)=x10;
x2(1)=x20;
y1(1)=y10;
y2(1)=y20;

% Zadana trajektoria zakłóceń
w1(1:n-1)=w10;
w2(1:n-1)=w20;

% Zadana trajektoria wyjść
y1zad(1:499)=y10;
y1zad(500:n)=65;
y2zad(1:999)=y20;
y2zad(1000:n)=25;

% Regulatory PID -----------------------(trzeba dostroić)-----------------
% classPID(K, Ti, Kd, Td, Tp, Hlim, Llim, Dir, AutoMan, ManVal) 
pid12=classPID(2, 10, 0, 0, 1, 100, 0, 1, 1, 0);
pid21=classPID(2, 10, 0, 0, 1, 100, 0, 1, 1, 0);

% Bloki odsprzęgające ----------------(trzeba sprawdzić)-----------------
% classLEADLAG(K, LEAD, LAG, Tp, Hlim, Llim) 
D21=classLEADLAG(-1/0.00007096, 0.02578/0.00007096, 0, 1, 100, 0);
D12=classLEADLAG(-1, 0, 0, 1, 100, 0);
du1=zeros(n-1,1); du2=zeros(n-1,1);

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
    
    if i<n
        % Sterowania wyznaczone przez regulatory PID
        uR2(i) =  pid12.calc(y1(i),y1zad(i));
        uR1(i) =  pid21.calc(y2(i),y2zad(i));
        
        % Człony odsprzęgające
        if czy_odsprz
            du2(i) = D21.calc(uR1(i));
            du1(i) = D12.calc(uR2(i));
        end
        
        % Sterowania po odsprzęganiu
        u2(i) = uR2(i) + du2(i);
        u1(i) = uR1(i) + du1(i);
        
        if u2(i)<0
            u2(i)=0;
        end
        if u1(i)<0
            u1(i)=0;
        end
    end

end    
 

% Wykresy
f=figure;
if czy_zlin 
    f.Name='Symulacja z modelem zlinearyzowanym';
else
    f.Name='Symulacja z modelem nieliniowym';
end
subplot(2,2,1)
stairs(tu, [u10; u1])
grid on;
title("Sygnal wejsciowy $u_1$",'interpreter','latex');
xlabel("Numer kroku",'interpreter','latex');
ylabel("$u_1$ $[cm^{3}/s]$",'interpreter','latex');

subplot(2,2,2)
stairs(tu, [u20; u2])
grid on;
title("Sygnal wejsciowy $u_2$",'interpreter','latex');
xlabel("Numer kroku",'interpreter','latex');
ylabel("$u_2$ $[cm^{3}/s]$",'interpreter','latex');

subplot(2,2,3);
plot(ty, [y1])
grid on;
hold on;
plot(ty, [y1zad]);
title("Sygnal wyjsciowy $y_1$",'interpreter','latex');
xlabel("Numer kroku",'interpreter','latex');
ylabel("$y_1, y_1^{zad}$ $[cm]$",'interpreter','latex');
% ylim([0.9*y1(1); 1.1*y1(1)]);

subplot(2,2,4);
plot(ty, [y2])
grid on;
hold on;
plot(ty, [y2zad]);
title("Sygnal wyjsciowy $y_2$",'interpreter','latex');
xlabel("Numer kroku",'interpreter','latex');
ylabel("$y_2, y_2^{zad}$ $[^\circ C]$",'interpreter','latex');
% ylim([0.9*y2(1); 1.1*y2(1)]);

