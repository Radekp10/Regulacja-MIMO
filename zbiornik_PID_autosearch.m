clear;

[punkt_pracy, parametry_modelu] = zbiornik_param;
[x10, x20, u10, u20, w10, w20, y10, y20] = punkt_pracy();
Yzad = zbiornik_y_zad(10000, 10, 20, [y10, y20]);

czy_zlin = false; % czy model liniowy czy zlinearyzowany
czy_odsprz = false; % czy PID z odsprzęganiem czy bez


fminfun = @(x) control_error(getY(x, Yzad, 0, czy_zlin, czy_odsprz), Yzad);
options = optimoptions('fmincon', 'Display', 'final-detailed');
[A, b] = PID_limits(2);
attemts = 3;

% Best settings searching
tic
for i = 1:attemts
    disp(i)
    stop(i).nastawy = fmincon(fminfun, rand(8, 1)*100, A, b, [], [], [], [], [], options);
    [Y, ~] = symulacja_PID_fun(Yzad, to_nastawy(stop(i).nastawy, 2), 0, czy_zlin, czy_odsprz);
    err(i) = control_error(Y, Yzad);
end
disp(datestr(datenum(0, 0, 0, 0, 0, toc), 'HH:MM:SS'));

% Show results
[~, I] = min(err);
regulators = to_nastawy(stop(I).nastawy, 2);
[Y, U] = symulacja_PID_fun(Yzad, regulators, 0, czy_zlin, czy_odsprz);

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

function Y = getY(x, Yzad, decoupling, czy_zlin, czy_odsprz)
[Y, ~] = symulacja_PID_fun(Yzad, to_nastawy(x, 2), decoupling, czy_zlin, czy_odsprz);
end

function reg_lokalne = to_nastawy(x, n)
for i = 1:n
    reg_lokalne(i).K = x((i - 1)*4+1);
    reg_lokalne(i).Ti = x((i - 1)*4+2);
    reg_lokalne(i).Kd = x((i - 1)*4+3);
    reg_lokalne(i).Td = x((i - 1)*4+4);
end
end

function [A, b] = PID_limits(n)
% Ax <= b
% x = [K; Ti; Td]
A = [];
b = [];
for i = 1:n
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

    K_b_up = 5000;
    K_b_down = 0.001;

    Ti_b_up = 5000;
    Ti_b_down = 0.001;

    Kd_b_up = 5000;
    Kd_b_down = 0.001;

    Td_b_up = 5000;
    Td_b_down = 0.001;

    pom = [K_A_up; K_A_down; Ti_A_up; Ti_A_down; Kd_A_up; Kd_A_down; Td_A_up; Td_A_down];
    A = [A; pom];

    pom = [K_b_up; K_b_down; Ti_b_up; Ti_b_down; Kd_b_up; Kd_b_down; Td_b_up; Td_b_down];
    b = [b; pom];
end
end
