[n, ~] = size(Yzad);
tu = 0:1:(n - 1);
ty = 1:1:n;

subplot(2, 2, 1)
stairs(tu, [u10; Unbi(:, 1)], 'Color', kolor(1))
hold on;
stairs(tu, [u10; Unbif(:, 1)], 'Color', kolor(3))
stairs(tu, [u10; Unb(:, 1)], 'Color', kolor(4))
stairs(tu, [u10; Unbf(:, 1)], 'Color', kolor(5))
% stairs(tu, [u10; Unzf(:, 1)], 'Color', kolor(3))
% stairs(tu, [u10; Ulbf(:, 1)], 'Color', kolor(4))
% stairs(tu, [u10; Ulzf(:, 1)], 'Color', kolor(5))
grid on;
title("Sygnal wejsciowy $u_1$", 'interpreter', 'latex');
xlabel("Numer kroku", 'interpreter', 'latex');
ylabel("$u_1$ $[cm^{3}/s]$", 'interpreter', 'latex');

subplot(2, 2, 2)
stairs(tu, [u20; Unbi(:, 2)], 'Color', kolor(1), 'DisplayName', 'Reg. PI, nastawy Z-N')
hold on;
stairs(tu, [u20; Unbif(:, 2)], 'Color', kolor(3), 'DisplayName', 'Reg. PI, nastawy fmincon')
stairs(tu, [u20; Unb(:, 2)], 'Color', kolor(4), 'DisplayName', 'Reg. PID, nastawy rêczne')
stairs(tu, [u20; Unbf(:, 2)], 'Color', kolor(5), 'DisplayName', 'Reg. PID, nastawy fmincon')
% stairs(tu, [u20; Unbf(:, 2)], 'Color', kolor(1), 'DisplayName', 'Model nieliniowy bez odsprzêgania')
% hold on;
% stairs(tu, [u20; Unzf(:, 2)], 'Color', kolor(3), 'DisplayName', 'Model nieliniowy z odprzêganiem')
% stairs(tu, [u20; Ulbf(:, 2)], 'Color', kolor(4), 'DisplayName', 'Model liniowy bez odsprzêgania')
% stairs(tu, [u20; Ulzf(:, 2)], 'Color', kolor(5), 'DisplayName', 'Model liniowy z odsprzêganiem')
grid on;
title("Sygnal wejsciowy $u_2$", 'interpreter', 'latex');
xlabel("Numer kroku", 'interpreter', 'latex');
ylabel("$u_2$ $[cm^{3}/s]$", 'interpreter', 'latex');
legend('show', 'location', 'best');

subplot(2, 2, 3);
plot(ty, Yzad(:, 1), 'Color', kolor(2));
hold on;
grid on;
stairs(ty, Ynbi(:, 1), 'Color', kolor(1))
stairs(ty, Ynbif(:, 1), 'Color', kolor(3))
stairs(ty, Ynb(:, 1), 'Color', kolor(4))
stairs(ty, Ynbf(:, 1), 'Color', kolor(5))
title("Sygnal wyjsciowy $y_1$", 'interpreter', 'latex');
xlabel("Numer kroku", 'interpreter', 'latex');
ylabel("$y_1, y_1^{zad}$ $[cm]$", 'interpreter', 'latex');
% ylim([0.9*y1(1); 1.1*y1(1)]);

subplot(2, 2, 4);
plot(ty, Yzad(:, 2), 'Color', kolor(2));
grid on;
hold on;
stairs(ty, Ynbi(:, 2), 'Color', kolor(1))
stairs(ty, Ynbif(:, 2), 'Color', kolor(3))
stairs(ty, Ynb(:, 2), 'Color', kolor(4))
stairs(ty, Ynbf(:, 2), 'Color', kolor(5))
title("Sygnal wyjsciowy $y_2$", 'interpreter', 'latex');
xlabel("Numer kroku", 'interpreter', 'latex');
ylabel("$y_2, y_2^{zad}$ $[^\circ C]$", 'interpreter', 'latex');
% ylim([0.9*y2(1); 1.1*y2(1)]);