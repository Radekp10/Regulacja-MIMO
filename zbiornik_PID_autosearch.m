clear;

[punkt_pracy, parametry_modelu] = zbiornik_param;
[x10, x20, u10, u20, w10, w20, y10, y20] = punkt_pracy();
pp.y1 = y10;
pp.y2 = y20;
pp.w1 = w10;
pp.w2 = w20;
% [Yzad ,W]= zbiornik_y_zad(10000, 10, 20, pp);
[Yzad, W] = zbiornik_y_zad(8000, 10, -1, pp);

czy_zlin = false; % czy model liniowy czy zlinearyzowany
czy_odsprz = false; % czy PID z odsprzęganiem czy bez
czy_rozniczka = true;


fminfun = @(x) sum(control_error(getY(x, Yzad, W, 0, czy_zlin, czy_odsprz, czy_rozniczka), Yzad));
options = optimoptions('fmincon', 'Display', 'off');
[A, b] = PID_limits(2, czy_rozniczka);
attemts = 3;

% Best settings searching
tic
for i = 1:attemts
    disp(i)
    if czy_rozniczka
        stop(i).nastawy = fmincon(fminfun, [1.5, 179, 44.75/4, 0, 1.44, 92, 23/8, 0], A, b, [], [], [], [], [], options);
    else
        stop(i).nastawy = fmincon(fminfun, [1.125, 300, 1.08, 150], A, b, [], [], [], [], [], options);
    end
    [Y, ~] = symulacja_PID_fun(Yzad, W, to_nastawy(stop(i).nastawy, 2, czy_rozniczka), 0, czy_zlin, czy_odsprz);
    err(i) = sum(control_error(Y, Yzad));
end
disp(datestr(datenum(0, 0, 0, 0, 0, toc), 'HH:MM:SS'));

% Show results
[~, I] = min(err);
regulators = to_nastawy(stop(I).nastawy, 2, czy_rozniczka);
[Y, U] = symulacja_PID_fun(Yzad, W, regulators, 0, czy_zlin, czy_odsprz);
control_error(Y, Yzad)

% Wykresy
f = figure;
if czy_zlin
    f.Name = 'Symulacja z modelem zlinearyzowanym';
else
    f.Name = 'Symulacja z modelem nieliniowym';
end

[n, ~] = size(Yzad);
tu = 0:1:(n - 1);
ty = 1:1:n;

subplot(2, 2, 1)
stairs(tu, [u10; U(:, 1)])
grid on;
title("Sygnal wejsciowy $u_1$", 'interpreter', 'latex');
xlabel("Numer kroku", 'interpreter', 'latex');
ylabel("$u_1$ $[cm^{3}/s]$", 'interpreter', 'latex');

subplot(2, 2, 2)
stairs(tu, [u20; U(:, 2)])
grid on;
title("Sygnal wejsciowy $u_2$", 'interpreter', 'latex');
xlabel("Numer kroku", 'interpreter', 'latex');
ylabel("$u_2$ $[cm^{3}/s]$", 'interpreter', 'latex');

subplot(2, 2, 3);
plot(ty, Y(:, 1))
grid on;
hold on;
plot(ty, Yzad(:, 1));
title("Sygnal wyjsciowy $y_1$", 'interpreter', 'latex');
xlabel("Numer kroku", 'interpreter', 'latex');
ylabel("$y_1, y_1^{zad}$ $[cm]$", 'interpreter', 'latex');
% ylim([0.9*y1(1); 1.1*y1(1)]);

subplot(2, 2, 4);
plot(ty, Y(:, 2))
grid on;
hold on;
plot(ty, Yzad(:, 2));
title("Sygnal wyjsciowy $y_2$", 'interpreter', 'latex');
xlabel("Numer kroku", 'interpreter', 'latex');
ylabel("$y_2, y_2^{zad}$ $[^\circ C]$", 'interpreter', 'latex');
% ylim([0.9*y2(1); 1.1*y2(1)]);

function Y = getY(x, Yzad, W, decoupling, czy_zlin, czy_odsprz, czy_rozniczka)
[Y, ~] = symulacja_PID_fun(Yzad, W, to_nastawy(x, 2, czy_rozniczka), decoupling, czy_zlin, czy_odsprz);
end

function reg_lokalne = to_nastawy(x, n, czy_rozniczka)
for i = 1:n
    if czy_rozniczka
        reg_lokalne(i).K = x((i - 1)*4+1);
        reg_lokalne(i).Ti = x((i - 1)*4+2);
        reg_lokalne(i).Kd = x((i - 1)*4+3);
        reg_lokalne(i).Td = x((i - 1)*4+4);
    else
        reg_lokalne(i).K = x((i - 1)*2+1);
        reg_lokalne(i).Ti = x((i - 1)*2+2);
        reg_lokalne(i).Kd = 0;
        reg_lokalne(i).Td = 0;
    end
end
end

function [A, b] = PID_limits(n, czy_rozniczka)
% Ax <= b
% x = [K; Ti; Td]
A = [];
b = [];
for i = 1:n
    if czy_rozniczka
        K_A_up = zeros(1, n*4);
        K_A_up((i - 1)*4+1) = 1;
        K_A_down = zeros(1, n*4);
        K_A_down((i - 1)*4+1) = -1;

        Ti_A_up = zeros(1, n*4);
        Ti_A_up((i - 1)*4+2) = 1;
        Ti_A_down = zeros(1, n*4);
        Ti_A_down((i - 1)*4+2) = -1;


        Kd_A_up = zeros(1, n*4);
        Kd_A_up((i - 1)*4+3) = 1;
        Kd_A_down = zeros(1, n*4);
        Kd_A_down((i - 1)*4+3) = -1;

        Td_A_up = zeros(1, n*4);
        Td_A_up((i - 1)*4+4) = 1;
        Td_A_down = zeros(1, n*4);
        Td_A_down((i - 1)*4+4) = -1;
    else
        K_A_up = zeros(1, n*2);
        K_A_up((i - 1)*2+1) = 1;
        K_A_down = zeros(1, n*2);
        K_A_down((i - 1)*2+1) = -1;

        Ti_A_up = zeros(1, n*2);
        Ti_A_up((i - 1)*2+2) = 1;
        Ti_A_down = zeros(1, n*2);
        Ti_A_down((i - 1)*2+2) = -1;
    end

    K_b_up = 5000;
    K_b_down = 0.001;

    Ti_b_up = 5000;
    Ti_b_down = 0.001;

    if czy_rozniczka
        Kd_b_up = 5000;
        Kd_b_down = 0.001;

        Td_b_up = 5000;
        Td_b_down = 0.001;
    end

    if czy_rozniczka
        pom = [K_A_up; K_A_down; Ti_A_up; Ti_A_down; Kd_A_up; Kd_A_down; Td_A_up; Td_A_down];
        A = [A; pom];

        pom = [K_b_up; K_b_down; Ti_b_up; Ti_b_down; Kd_b_up; Kd_b_down; Td_b_up; Td_b_down];
        b = [b; pom];
    else
        pom = [K_A_up; K_A_down; Ti_A_up; Ti_A_down];
        A = [A; pom];

        pom = [K_b_up; K_b_down; Ti_b_up; Ti_b_down];
        b = [b; pom];
    end
end
end
