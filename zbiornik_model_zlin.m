% Zlinearyzowany model zbiornika

function [f1Lh, f2Lh, gLh] = zbiornik_model_zlin
    f1Lh=@f1L;
    f2Lh=@f2L;
    gLh=@gL;
end


% Zbiornik z mieszaniem  - model w przestrzeni stanu: dynamiczny, ciągły i zlinearyzowany
function [f1] = f1L(u1, u2, w1, w2, x1, x2)
    [punkt_pracy, parametry_modelu] = zbiornik_param;
    [T_H,T_C,tau,tau_C,C,alpha] = parametry_modelu();
    [x10, x20, u10, u20]=punkt_pracy();
    % Równanie zlinearyzowane:
%     f1=1/(2*C)*(u10/x10-u10/(x10^2)*(x1-x10)+1/x10*(u1-u10)) + 1/(2*C)*(u20/x10-u20/(x10^2)*(x1-x10)+1/x10*(u2-u20))  +  w1/(2*C)*((1/x10)-(1/x10^2)*(x1-x10)) -alpha/(2*C)*((1/sqrt(x10))-(1/(2*sqrt(x10^3)))*(x1-x10));
    
    % Równanie zlinearyzowane (to samo co wyżej) po uporządkowaniu
    % zmiennych do postaci: f1=a11*x1 + 0*x2 + b11*u1 + b12*u2 + const
    f1=1/(2*C)*(-u10/(x10^2)-u20/(x10^2)-w1/(x10^2)-alpha/(2*sqrt(x10^3)))*x1 + 1/(2*C)*(1/x10)*u1 + 1/(2*C)*(1/x10)*u2 + 1/(2*C)*((u10+u20+2*w1)/x10 - alpha/(2*sqrt(x10)));
end

function [f2] = f2L(u1, u2, w1, w2, x1, x2)
    [punkt_pracy, parametry_modelu] = zbiornik_param;
    [T_H,T_C,tau,tau_C,C,alpha] = parametry_modelu();
    [x10, x20, u10, u20]=punkt_pracy();
    % Równanie zlinearyzowane:
%     f2=T_H/C*(u10/(x10^2)-2*u10/(x10^3)*(x1-x10)+1/(x10^2)*(u1-u10))  +  T_C/C*(u20/(x10^2)-2*u20/(x10^3)*(x1-x10)+1/(x10^2)*(u2-u20))  +  (w1*w2)/C*((1/x10^2)-(2/x10^3)*(x1-x10)) -alpha/C*(x20/sqrt(x10^3)+(1/sqrt(x10^3))*(x2-x20)-(3*x20/(2*sqrt(x10^5)))*(x1-x10));
    
    % Równanie zlinearyzowane (to samo co wyżej) po uporządkowaniu
    % zmiennych do postaci: f2=a21*x1 + a22*x2 + b21*u1 + b22*u2 + const
    f2=1/C*(-2*T_H*u10/(x10^3)-2*T_C*u20/(x10^3)-2*w2*w1/(x10^3)-alpha*(-3)*x20/(2*sqrt(x10^5)))*x1 + -alpha/C*(1/sqrt(x10^3))*x2 + T_H/C*(1/(x10^2))*u1 + T_C/C*(1/(x10^2))*u2 + T_H/C*2*u10/(x10^2) + T_C/C*2*u20/(x10^2) + w1*w2*3/(C*x10^2) - alpha*3*x20/(C*2*sqrt(x10^3));
end

function [y] = gL(x1, x2, x2_opoz)
    y(1)=x1;
    y(2)=x2_opoz;
end   
