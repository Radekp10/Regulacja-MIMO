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
    [x10, x20]=punkt_pracy();
    f1=(u1+u2+w1)/(2*C)*((1/x10)-(1/x10^2)*(x1-x10)) -alpha/(2*C)*((1/sqrt(x10))-(1/(2*sqrt(x10^3)))*(x1-x10));
end

function [f2] = f2L(u1, u2, w1, w2, x1, x2)
    [punkt_pracy, parametry_modelu] = zbiornik_param;
    [T_H,T_C,tau,tau_C,C,alpha] = parametry_modelu();
    [x10, x20]=punkt_pracy();
    f2=(u1*T_H + u2*T_C + w1*w2)/C*((1/x10^2)-(2/x10^3)*(x1-x10)) -alpha/C*(x20/sqrt(x10^3)+(1/sqrt(x10^3))*(x2-x20)-(3*x20/(2*sqrt(x10^5)))*(x1-x10));
end

function [y] = gL(x1, x2, x2_opoz)
    y(1)=x1;
    y(2)=x2_opoz;
end   
