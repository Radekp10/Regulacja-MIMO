% Skrypt do regulacji zbiornikiem przy pomocy analitycznego DMC 
% z przycinaniem sterowania

% clear;

% Pobranie uchwytu do funkcji symulującej zbiornik
[punkt_pracy, parametry_modelu]=zbiornik_param;
[symulacja_obiektu]=zbiornik_sym;

% Nominalny punkt pracy
[x10, x20, u10, u20, w10, w20, y10, y20]=punkt_pracy();

% Parametry modelu
[T_H,T_C,tau,tau_C,C,alpha] = parametry_modelu();

Tp=10; % okres próbkowania obiektu (i działania regulatora) [s]
h=1; % wielkość kroku (czyli jak szybko chodzi sam obiekt) [s]
n=8e3/Tp; % liczba kroków 'Tp'

czy_zlin=false; % czy model nieliniowy czy zlinearyzowany
czy_komp_zakl=true; % czy kompensować zakl. mierzone

% Osie do wykresu
tu=0:Tp:Tp*(n-1);
ty=Tp:Tp:Tp*n;

x1=zeros(n*Tp/h,1);
x2=zeros(n*Tp/h,1);
y1=zeros(n,1);
y2=zeros(n,1);
y1zad=zeros(n,1);
y2zad=zeros(n,1);
e1=0; % uchyby regulacji
e2=0;
u1=zeros(n-1,1);
u2=zeros(n-1,1);
w1=zeros(n-1,1);
w2=zeros(n-1,1);

% Początkowy punkt pracy
x1(1:Tp/h)=x10;
x2(1:Tp/h)=x20;
y1(1)=y10;
y2(1)=y20;

% Zadana trajektoria zakłóceń
w1(1:n-1)=w10;
w1(5000/Tp:n-1)=15;
w2(1:n-1)=w20;
w2(7000/Tp:n-1)=25;

% Zadana trajektoria wyjść
y1zad(1:floor(999/Tp))=y10;
y1zad(1000/Tp:n)=60;
y2zad(1:floor(2999/Tp))=y20;
y2zad(3000/Tp:n)=30;

% Horyzonty
D=200;
Dz=200;
N=30;
Nu=15;

% Macierzowa odp. skokowa
S=zbiornik_mac_odp_skok(Tp,Tp*D);
[ny,nu,~]=size(S);
Sz=zbiornik_mac_odp_skok_zakl(Tp,Tp*Dz);
[~,nz,~]=size(Sz);

% Wspolczynniki 
lambda1=1; 
lambda2=1;
psi1=1; psi2=1;
lambda=diag([lambda1, lambda2]);
psi=diag([psi1, psi2]);
Lambda=zeros(nu*Nu,nu*Nu);
Psi=zeros(ny*N,ny*N);
for j=1:Nu
    Lambda(1+(j-1)*nu:j*nu, 1+(j-1)*nu:j*nu)=lambda;
end
for j=1:N
    Psi(1+(j-1)*ny:j*ny, 1+(j-1)*ny:j*ny)=psi;
end


% Ograniczenia sterowania
u1min=0;
u1max=20;
u2min=30;
u2max=50;
u1max_ogr=ones(n,1)*u1max; % sygnal do wykresu
u1min_ogr=ones(n,1)*u1min; % sygnal do wykresu
u2max_ogr=ones(n,1)*u2max; % sygnal do wykresu
u2min_ogr=ones(n,1)*u2min; % sygnal do wykresu
u1min=u1min-u10;
u1max=u1max-u10;
u2min=u2min-u20;
u2max=u2max-u20;
du1max=10;
du2max=10;

% Ograniczenia wyjścia
y1min=0;
y1max=100;  
y2min=0;
y2max=100; 
y1min_ogr=ones(n,1)*y1min; % sygnal do wykresu
y1max_ogr=ones(n,1)*y1max; % sygnal do wykresu
y2min_ogr=ones(n,1)*y2min; % sygnal do wykresu
y2max_ogr=ones(n,1)*y2max; % sygnal do wykresu


% Wektor przeszłych przyrostow sterowania i zakłóceń
deltaupk=zeros((D-1)*nu,1);
deltazpk=zeros(Dz*nz,1);

% Generacja macierzy M i MP i MzP
[M,MP]=macierzeDMC(S,N,Nu);
MzP=macierzMzP(Sz,N);  % (N*ny) x (Dz*nz)

% Obliczanie parametrów regulatora
K=((M'*Psi*M+Lambda)^-1)*M'*Psi;
K1=K(1:nu,:);  % (nu) x (ny*N)
Ke=zeros(nu,ny);
for p=1:N
    K1p(p)={K1(:,1+(p-1)*ny:p*ny)};  % (nu) x (ny)
    Ke=Ke+K1p{p};  % (nu) x (ny)
end
for j=1:D-1
    Kuj(j)={K1*MP(:,1+(j-1)*nu:j*nu)};  % (nu) x (nu)
end
for j=1:Dz
    Kzj(j)={K1*MzP(:,1+(j-1)*nz:j*nz)};  % (nu) x (nz)
end


for i=2:n % czas symulacji

    if i-1-tau_C/Tp > 0
        u2_opoz=u2(i-1-tau_C/Tp);
    else
        u2_opoz=u2(1);
    end
    
    % Symulacja obiektu
    for j=1:Tp/h
        if (i-1)*Tp/h+j-1-tau/h > 0
            x2_opoz=x2((i-1)*Tp/h+j-1-tau/h);
        else
            x2_opoz=x2(1);
        end
        [y, x_nast] = symulacja_obiektu(u1(i-1)+u10, u2_opoz+u20, w1(i-1), w2(i-1), x1((i-1)*Tp/h+j-1), x2((i-1)*Tp/h+j-1), x2_opoz, h, czy_zlin);
        x1((i-1)*Tp/h+j)=x_nast(1);
        x2((i-1)*Tp/h+j)=x_nast(2);
    end
    y1(i)=y(1);
    y2(i)=y(2);
    
        
    % Symulację kończymy po wykonaniu ostatniego pomiaru wyjścia y(kk)
    % (u(kk) nie jest liczone, bo jego efekt byłby widoczny 
    % dopiero w y(kk+1) czyli już po końcu symulacji)
    if i==n
        break;
    end    
        
        
    for j=Dz:-1:2
       % przesuniecie wartosci w wektorze przeszlych przyrostow zakl.
       deltazpk(1+(j-1)*nz:j*nz)=deltazpk(1+(j-2)*nz:(j-1)*nz);    
    end
    deltazpk(1)=w1(i)-w1(i-1);
    deltazpk(2)=w2(i)-w2(i-1);         

    % Wyjście w chwili k
    Yk=[y1(i); y2(i)];       

    % Trajektoria refencyjna stała na horyzoncie predykcji
    Yref=[y1zad(i); y2zad(i)];

%        % Wersja alternatywna
%        Kuj_sum=zeros(nu,1);
%        for j=1:D-1
%            Kuj_sum=Kuj_sum+Kuj{j}*deltaupk(1+(j-1)*nu:j*nu);
%        end
%        % Sprawdzenie poprawności
%        if abs(Kuj_sum - K1*MP*deltaupk)>1e-6*ones(2,1), disp('err'); end

    if czy_komp_zakl==1
        deltauk=Ke*(Yref-Yk)-K1*MP*deltaupk-K1*MzP*deltazpk;
    elseif czy_komp_zakl==0
        deltauk=Ke*(Yref-Yk)-K1*MP*deltaupk;
    end

    if deltauk(1)>du1max
        deltauk(1)=du1max;
    end
    if deltauk(1)<-du1max
        deltauk(1)=-du1max;
    end    

    if deltauk(2)>du2max
        deltauk(2)=du2max;
    end
    if deltauk(2)<-du2max
        deltauk(2)=-du2max;
    end  
       
    for j=D-1:-1:2
       % przesuniecie wartosci w wektorze przeszlych przyrostow sterowania
       deltaupk(1+(j-1)*nu:j*nu)=deltaupk(1+(j-2)*nu:(j-1)*nu);    
    end
    deltaupk(1:2)=deltauk;

    u1(i)=u1(i-1)+deltauk(1); % nowe sterowania
    u2(i)=u2(i-1)+deltauk(2);

    if u1(i)>u1max
        u1(i)=u1max;
    end
    if u1(i)<u1min
        u1(i)=u1min;
    end

    if u2(i)>u2max
        u2(i)=u2max;
    end
    if u2(i)<u2min
        u2(i)=u2min;
    end
    deltaupk(1)=u1(i)-u1(i-1); % jesli u(k) zahaczylo o jakies ograniczenie, to trzeba jeszcze raz policzyc deltaupk(1)
    deltaupk(2)=u2(i)-u2(i-1);

  % Wskaźniki jakości regulacji
    e1 = e1 + (y1zad(i)-y1(i))^2;
    e2 = e2 + (y2zad(i)-y2(i))^2;
        
end
    
fprintf('E1=%.2f\n',e1');
fprintf('E2=%.2f\n',e2');

% Wykresy
close all;
f=figure;
if czy_zlin 
    f.Name='Symulacja z modelem zlinearyzowanym';
else
    f.Name='Symulacja z modelem nieliniowym';
end
subplot(2,2,1)
stairs(tu, [u10; u1+u10])
grid on;
hold on;
stairs(tu, u1max_ogr, '-.r');
stairs(tu, u1min_ogr, '-.k');
title("Sygnal wejsciowy $u_1$",'interpreter','latex');
xlabel("Czas [s]",'interpreter','latex');
ylabel("$u_1$ $[cm^{3}/s]$",'interpreter','latex');
legend('$u_1$', '$u^\mathrm{max}_1$', '$u^\mathrm{min}_1$','interpreter','latex','location','best');
ylim([u1min+u10-5; u1max+u10+5]);

subplot(2,2,2)
stairs(tu, [u20; u2+u20])
grid on;
hold on;
stairs(tu, u2max_ogr, '-.r');
stairs(tu, u2min_ogr, '-.k');
title("Sygnal wejsciowy $u_2$",'interpreter','latex');
xlabel("Czas [s]",'interpreter','latex');
ylabel("$u_2$ $[cm^{3}/s]$",'interpreter','latex');
legend('$u_2$', '$u^\mathrm{max}_2$', '$u^\mathrm{min}_2$','interpreter','latex','location','best');
ylim([u2min+u20-5; u2max+u20+5]);

subplot(2,2,3);
plot(ty, [y1])
grid on;
hold on;
plot(ty, [y1zad]);
title("Sygnal wyjsciowy $y_1$",'interpreter','latex');
xlabel("Czas [s]",'interpreter','latex');
ylabel("$y_1, y_1^{zad}$ $[cm]$",'interpreter','latex');
legend('$y_1$', '$y^\mathrm{zad}_1$','interpreter','latex','location','best'); 
% ylim([0.9*y1(1); 1.1*y1(1)]);

subplot(2,2,4);
plot(ty, [y2])
grid on;
hold on;
plot(ty, [y2zad]);
title("Sygnal wyjsciowy $y_2$",'interpreter','latex');
xlabel("Czas [s]",'interpreter','latex');
ylabel("$y_2, y_2^{zad}$ $[^\circ C]$",'interpreter','latex');
legend('$y_2$', '$y^\mathrm{zad}_2$','interpreter','latex','location','best'); 
% ylim([0.9*y2(1); 1.1*y2(1)]);


