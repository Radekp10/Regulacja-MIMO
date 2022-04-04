% Skrypt do symulacji i porównania modeli zbiornika
% Skrypt pozwala zmieniaæ wartoœci zmiennych steruj¹cych (u1 lub u2) i po
% restarcie nanosi kolejne linie wykres

clear;

ktore_sterowanie = 1; % czy zmieniamy sterowanie u1 czy u2

if ktore_sterowanie == 1
    u_list = [7, 9, 11, 13, 15, 17, 19];
elseif ktore_sterowanie == 2
    u_list = [20, 25, 30, 35, 40, 45, 50];
else
    ktore_sterowanie = 1;
    u_list = 13;
end

% Pobranie uchwytu do funkcji symulujÄ…cej zbiornik
[punkt_pracy, parametry_modelu] = zbiornik_param;
[symulacja_obiektu] = zbiornik_sym;

% Nominalny punkt pracy
[x10, x20, u10, u20, w10, w20, y10, y20] = punkt_pracy();

% Parametry modelu
[T_H, T_C, tau, tau_C, C, alpha] = parametry_modelu();

n = 2e3; % liczba krokÃ³w h
h = 1; % wielkoÅ›Ä‡ kroku [s]

% Osie do wykresu
tu = 0:1:(n - 1);
ty = 1:1:n;

x1 = zeros(n, 1);
x1l = zeros(n, 1);
x2 = zeros(n, 1);
x2l = zeros(n, 1);
y1 = zeros(n, 1);
y1l = zeros(n, 1);
y2 = zeros(n, 1);
y2l = zeros(n, 1);
u1 = zeros(n-1, 1);
u2 = zeros(n-1, 1);
w1 = zeros(n-1, 1);
w2 = zeros(n-1, 1);

% PoczÄ…tkowy punkt pracy
x1(1) = x10;
x1l(1) = x10;
x2(1) = x20;
x2l(1) = x20;
y1(1) = y10;
y1l(1) = y10;
y2(1) = y20;
y2l(1) = y20;


for i_list = 1:length(u_list)

    % Zadana trajektoria sterowaÅ„ i zakÅ‚Ã³ceÅ„
    if ktore_sterowanie == 1
        u1(1) = u10;
        u1(2:n-1) = u_list(i_list);
        u2(1:n-1) = u20;
    else
        u1(1:n-1) = u10;
        u2(1) = u20;
        u2(2:n-1) = u_list(i_list);
    end
    w1(1:n-1) = w10;
    w2(1:n-1) = w20;


    % GÅ‚Ã³wna pÄ™tla symulacji obiektu
    for i = 2:n
        if i - 1 - tau_C > 0
            u2_opoz = u2(i-1-tau_C);
        else
            u2_opoz = u2(1);
        end

        if i - 1 - tau > 0
            x2_opoz = x2(i-1-tau);
            x2_opozl = x2l(i-1-tau);
        else
            x2_opoz = x2(1);
            x2_opozl = x2l(1);
        end

        [y, x_nast] = symulacja_obiektu(u1(i-1), u2_opoz, w1(i-1), w2(i-1), x1(i-1), x2(i-1), x2_opoz, h, false);
        x1(i) = x_nast(1);
        x2(i) = x_nast(2);
        y1(i) = y(1);
        y2(i) = y(2);


        [yl, x_nastl] = symulacja_obiektu(u1(i-1), u2_opoz, w1(i-1), w2(i-1), x1l(i-1), x2l(i-1), x2_opozl, h, true);
        x1l(i) = x_nastl(1);
        x2l(i) = x_nastl(2);
        y1l(i) = yl(1);
        y2l(i) = yl(2);
    end


    % Wykresy
    f = figure(17); % Losowy numer, ale pozwala zawsze odnosiæ siê do tego samego wykresu
    f.Name = 'Porównanie symulacji modelu nieliniowego i zlinearyzowanego';

    if ktore_sterowanie == 1
        sterowanie = "u_1";
    elseif ktore_sterowanie == 2
        sterowanie = "u_2";
    else
        sterowanie = "u";
    end

    subplot(2, 1, 1)
    hold on;
    stairs(ty, y1, 'color', kolor(i_list), 'DisplayName', sprintf("%s = %d[cm^{3}/s]", sterowanie, u_list(i_list)))
    stairs(ty, y1l, '--', 'color', kolor(i_list), 'HandleVisibility', 'off')
    grid on;
    title(sprintf("Porownanie odpowiedzi $y_1$ na skok sterowania $%s$", sterowanie), 'interpreter', 'latex');
    xlabel("Numer kroku", 'interpreter', 'latex');
    ylabel("$y_1$ $[cm]$", 'interpreter', 'latex');
    legend('show', 'location', 'best');

    subplot(2, 1, 2)
    hold on;
    stairs(ty, y2, 'color', kolor(i_list))
    stairs(ty, y2l, '--', 'color', kolor(i_list))
    grid on;
    title(sprintf("Porownanie odpowiedzi $y_2$ na skok sterowania $%s$", sterowanie), 'interpreter', 'latex');
    xlabel("Numer kroku", 'interpreter', 'latex');
    ylabel("$y_2$ $[^\circ C]$", 'interpreter', 'latex');
    
%     [y1(2000), y1l(2000), abs(y1l(2000)-y1(2000)), 100*abs(y1l(2000)-y1(2000))/y1(2000)]
%     [y2(2000), y2l(2000), abs(y2l(2000)-y2(2000)), 100*abs(y2l(2000)-y2(2000))/y2(2000)]
%     pause;
end
