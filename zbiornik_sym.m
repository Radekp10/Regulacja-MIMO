% Funkcja do symulacji zbiornika z mieszaniem

function [symh] = zbiornik_sym
    symh=@symulacja_obiektu;
end

% Symulacja zbiornika - rozwiązanie numeryczne układu równań metodą
% Rungego-Kutty IV rzędu (metoda jawna, jednokrokowa) z krokiem 'h'
function [y, x_nast] = symulacja_obiektu(u1, u2, w1, w2, x1, x2, x2_opoz, h, czy_zlin)
    % u=[F_H, F_Cin];
    % x=[h, T];
    % w=[F_D, T_D];
    % y=[h, T_out];
    
    if czy_zlin
        [f1, f2, g] = zbiornik_model_zlin2;
    else
        [f1, f2, g] = zbiornik_model;
    end

    % K{i}{j} - współczynnik K{i} dla zmiennej x{j}

    K11=f1(u1, u2, w1, w2, x1, x2);
    K12=f2(u1, u2, w1, w2, x1, x2);

    K21=f1(u1, u2, w1, w2, x1+0.5*K11*h, x2+0.5*K12*h);
    K22=f2(u1, u2, w1, w2, x1+0.5*K11*h, x2+0.5*K12*h);

    K31=f1(u1, u2, w1, w2, x1+0.5*K21*h, x2+0.5*K22*h);
    K32=f2(u1, u2, w1, w2, x1+0.5*K21*h, x2+0.5*K22*h);

    K41=f1(u1, u2, w1, w2, x1+K31*h, x2+K32*h);
    K42=f2(u1, u2, w1, w2, x1+K31*h, x2+K32*h);

    x_nast(1)=x1+h/6*(K11+2*K21+2*K31+K41);
    x_nast(2)=x2+h/6*(K12+2*K22+2*K32+K42);

    y=g(x_nast(1), x_nast(2), x2_opoz);
    
end
