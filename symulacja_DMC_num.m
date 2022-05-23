% Skrypt do regulacji zbiornikiem przy pomocy numerycznego DMC

clear;

% Pobranie uchwytu do funkcji symulującej zbiornik
[punkt_pracy, parametry_modelu]=zbiornik_param;
[symulacja_obiektu]=zbiornik_sym;

% Nominalny punkt pracy
[x10, x20, u10, u20, w10, w20, y10, y20]=punkt_pracy();

% Parametry modelu
[T_H,T_C,tau,tau_C,C,alpha] = parametry_modelu();

% Uwaga: 'Tp' musi być wielokrotnością 'h'
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

Ymax=zeros(N*ny,1);
for j=1:N
   Ymax(1+(j-1)*ny:j*ny)=[y1max-y10, y2max-y20];
end

Ymin=zeros(N*ny,1);
for j=1:N
   Ymin(1+(j-1)*ny:j*ny)=[y1min-y10, y2min-y20];
end

Umax=zeros(Nu*nu,1);
for j=1:Nu
    Umax(1+(j-1)*nu:j*nu)=[u1max-u10, u2max-u20];
end

Umin=zeros(Nu*nu,1);
for j=1:Nu
    Umin(1+(j-1)*nu:j*nu)=[u1min-u10, u2min-u20];
end

dUmax=zeros(Nu*nu,1);
for j=1:Nu
    dUmax(1+(j-1)*nu:j*nu)=[du1max, du2max];
end

% Wektor przeszłych przyrostow sterowania i zakłóceń
deltaupk=zeros((D-1)*nu,1);
deltazpk=zeros(Dz*nz,1);

% Generacja macierzy M i MP i MzP
[M,MP]=macierzeDMC(S,N,Nu);
MzP=macierzMzP(Sz,N);  % (N*ny) x (Dz*nz)
  

% Głowna petla regulacji
for i=2:n

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
    
    if mod(i,100)==0
        fprintf('Wykonane iteracje: %d/%d\n',i,n);
    end
    
    for j=Dz:-1:2
          % przesuniecie wartosci w wektorze przeszlych przyrostow zakl.
          deltazpk(1+(j-1)*nz:j*nz)=deltazpk(1+(j-2)*nz:(j-1)*nz);    
    end
    deltazpk(1)=w1(i)-w1(i-1);
    deltazpk(2)=w2(i)-w2(i-1);


    Yk=zeros(N*ny,1);
    for j=1:N
        Yk(1+(j-1)*ny:j*ny)=[y1(i)-y10, y2(i)-y20];
    end

    Yzad=zeros(N*ny,1);
    for j=1:N
        Yzad(1+(j-1)*ny:j*ny)=[y1zad(i)-y10, y2zad(i)-y20];
    end

    J=tril(ones(Nu*nu)); % macierz trojkatna dolna
 
    Uk1=zeros(Nu*nu,1); % U(k-1)
    for j=1:Nu
        Uk1(1+(j-1)*nu:j*nu)=[u1(i-1), u2(i-1)];
    end

    if czy_komp_zakl==1
        Y0=Yk+MP*deltaupk+MzP*deltazpk;
    elseif czy_komp_zakl==0
        Y0=Yk+MP*deltaupk;
    end

    H=2*(M'*M+Lambda);
    f=-2*M'*(Yzad-Y0);
    A=[-J; J; -M; M];
    b=[-Umin+Uk1; Umax-Uk1; -Ymin+Y0; Ymax-Y0];

    options=optimoptions(@quadprog,'MaxIterations',200 ,'ConstraintTolerance',1e-8,'Display','off'); %default: 200, 1e-8
    [deltau,fval,exitflag]=quadprog(H,f,A,b,[],[],-dUmax,dUmax,[],options);

    if exitflag > 0
        for j=D-1:-1:2
            % przesuniecie wartosci w wektorze przeszlych przyrostow sterowania
            deltaupk(1+(j-1)*ny:j*ny)=deltaupk(1+(j-2)*ny:(j-1)*ny);    
        end
        deltaupk(1:2)=deltau(1:2);
        u1(i)=u1(i-1)+deltau(1); % nowe sterowania
        u2(i)=u2(i-1)+deltau(2);
    else
        fprintf('Exitflag for QP: %d in iter %d\n',exitflag,i);
        u1(i)=u1(i-1); % przyjęcie sterowania wyznaczonego w poprzedniej chwili probkowania
        u2(i)=u2(i-1);
    end

    % Wskaźniki jakości regulacji
    e1 = e1 + (y1zad(i)-y1(i))^2;
    e2 = e2 + (y2zad(i)-y2(i))^2;
       
end    


fprintf('E1=%.2f\n',e1');
fprintf('E2=%.2f\n',e2');

% Wykresy
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
ylim([u1min-5; u1max+5]);

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
ylim([u2min-5; u2max+5]);

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

