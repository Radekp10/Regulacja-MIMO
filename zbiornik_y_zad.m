function Yzad = zbiornik_y_zad(n, step_amount, seed, work_point)
% Funkcja wyznaczajaca zadana trajektorie zbiornika
% n - dlugosc wektora
% step_amount - ilosc skokow wartosci zadanych (losowo wybierane jest czy
% jedna, czy dwie wartosci sie zmieniaja
% seed - ziarno losowe
% work_point - punkt pracy, wektor 1x2: 1 - h, 2 - T

% Arbitralnie wybrane graniczne wartosci temperatry i wysokosci cieczy
h_max = 81;
h_min = 27;
T_max = 51;
T_min = 17;


Yzad = zeros(n, 2);
generator = rng(seed);
range = round(n/step_amount);
prev_work_point = work_point;

for i = 1:step_amount
    % Ustalenie ile wartosci ma skok
    switch randsample(3, 1)
        case 1
            Yzad((i - 1)*range+1:i*range, 1) = h_min + (h_max - h_min) * rand;
            Yzad((i - 1)*range+1:i*range, 2) = prev_work_point(2);
            prev_work_point(1) = Yzad(i * range,1);
        case 2
            Yzad((i - 1)*range+1:i*range, 1) = prev_work_point(1);
            Yzad((i - 1)*range+1:i*range, 2) = T_min + (T_max - T_min) * rand;
            prev_work_point(2) = Yzad(i * range,2);
        otherwise
            Yzad((i - 1)*range+1:i*range, 1) = h_min + (h_max - h_min) * rand;
            Yzad((i - 1)*range+1:i*range, 2) = T_min + (T_max - T_min) * rand;
            prev_work_point(1) = Yzad(i * range,1);
            prev_work_point(2) = Yzad(i * range,2);
    end
end
end