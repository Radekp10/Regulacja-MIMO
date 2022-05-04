function [Y, U] = symulacja_PID_fun(Yzad, W, regulators, decoupling, czy_zlin, czy_odsprz)


% Pobranie uchwytu do funkcji symulującej zbiornik
[punkt_pracy, parametry_modelu] = zbiornik_param;
[symulacja_obiektu] = zbiornik_sym;

% Nominalny punkt pracy
[x10, x20, u10, u20, w10, w20, y10, y20] = punkt_pracy();

% Parametry modelu
[T_H, T_C, tau, tau_C, C, alpha] = parametry_modelu();

[n, ~] = size(Yzad); % liczba kroków h
h = 1; % wielkość kroku [s]

x1 = zeros(n, 1);
x2 = zeros(n, 1);
y1 = zeros(n, 1);
y2 = zeros(n, 1);
u1 = zeros(n-1, 1);
u2 = zeros(n-1, 1);
uR1 = zeros(n-1, 1);
uR2 = zeros(n-1, 1);

% Początkowy punkt pracy
x1(1) = x10;
x2(1) = x20;
y1(1) = y10;
y2(1) = y20;
u1(1) = u10;
u2(1) = u20;

% Zadana trajektoria zakłóceń
w1 = W(1:n-1, 1);
w2 = W(1:n-1, 2);

% Zadana trajektoria wyjść
y1zad = Yzad(:, 1);
y1zad(1) = y10;
y2zad = Yzad(:, 2);
y2zad(1) = y20;

% Regulatory PID -----------------------(trzeba dostroić)-----------------
% classPID(K, Ti, Kd, Td, Tp, Hlim, Llim, Dir, AutoMan, ManVal)
pid12 = ovation_classPID(regulators(1).K, regulators(1).Ti, regulators(1).Kd, regulators(1).Td, 1, 100, 0, 1, 1, 0);
pid21 = ovation_classPID(regulators(2).K, regulators(2).Ti, regulators(2).Kd, regulators(2).Td, 1, 100, 0, 1, 1, 0);

% Bloki odsprzęgające
% classLEADLAG(K, LEAD, LAG, Tp, Hlim, Llim)
D21 = ovation_classLEADLAG(0.00002384/0.00007096, -0.01056/0.00002384, 0.02578/0.00007096, 1, 100, -100);
D12 = ovation_classLEADLAG(0, 0, 0, 1, 100, -100); %wyłączony
du1 = zeros(n-1, 1);
du2 = zeros(n-1, 1);

% Sterowanie ręczne na pocz. symulacji (tylko pierwsza iteracja głównej
% pętli, by poprawnie zainicjalizował się PID)
% SetAutoMan(obj, AutoMan, ManVal)
pid12.SetAutoMan(0, u20);
pid21.SetAutoMan(0, u10);


% Główna pętla symulacji obiektu
for i = 2:n
    if i - 1 - tau_C > 0
        u2_opoz = u2(i-1-tau_C);
    else
        u2_opoz = u2(1);
    end

    if i - 1 - tau > 0
        x2_opoz = x2(i-1-tau);
    else
        x2_opoz = x2(1);
    end

    [y, x_nast] = symulacja_obiektu(u1(i-1), u2_opoz, w1(i-1), w2(i-1), x1(i-1), x2(i-1), x2_opoz, h, czy_zlin);
    x1(i) = x_nast(1);
    x2(i) = x_nast(2);
    y1(i) = y(1);
    y2(i) = y(2);

    % Przejście na ster. automatyczne
    if i > 2
        pid12.SetAutoMan(1, u20);
        pid21.SetAutoMan(1, u10);
    end

    if i < n
        % Sterowania wyznaczone przez regulatory PID
        uR2(i) = pid12.calc(y1(i), y1zad(i));
        uR1(i) = pid21.calc(y2(i), y2zad(i));

        % Człony odsprzęgające
        if czy_odsprz
            du2(i) = D21.calc(uR1(i));
            du1(i) = D12.calc(uR2(i));
        end

        % Sterowania po odsprzęganiu
        u2(i) = uR2(i) + du2(i);
        u1(i) = uR1(i) + du1(i);

        if u2(i) < 0
            u2(i) = 0;
        end
        if u1(i) < 0
            u1(i) = 0;
        end

    end

end

Y = [y1, y2];
U = [u1, u2];