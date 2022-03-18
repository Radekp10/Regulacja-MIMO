% Funkcja zwracająca uchwyty do funkcji lokalnych
function [pph, pmh] = zbiornik_param
    pph=@punkt_pracy;
    pmh=@parametry_modelu;
end

% Nominalny punkt pracy zbiornika
function [x10, x20, u10, u20, w10, w20, y10, y20]=punkt_pracy
    F_C=35;
    F_H=13;
    F_D=11;
    T_D=31;
    h=54.39;
    T=34.05;
    
    x10=h;
    x20=T;
    u10=F_H;
    u20=F_C;
    w10=F_D;
    w20=T_D;
    y10=h;
    y20=T;
end

% Dodatkowe parametry modelu
function [T_H,T_C,tau,tau_C,C,alpha] = parametry_modelu
    C=0.8;
    alpha=8;
    T_C=25;
    T_H=61;
    tau_C=120;
    tau=60;
end

