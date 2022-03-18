% Dynamiczny, ciągły i nieliniowy model zbiornika

function [f1h, f2h, gh] = zbiornik_model
    f1h=@f1;
    f2h=@f2;
    gh=@g;
end


% Zbiornik z mieszaniem  - model w przestrzeni stanu
function [f1] = f1(u1, u2, w1, w2, x1, x2)
    [~, parametry_modelu] = zbiornik_param;
    [T_H,T_C,tau,tau_C,C,alpha] = parametry_modelu();
    f1=(u1+u2+w1-alpha*sqrt(x1)) / (2*C*x1);
end

function [f2] = f2(u1, u2, w1, w2, x1, x2)
    [~, parametry_modelu] = zbiornik_param;
    [T_H,T_C,tau,tau_C,C,alpha] = parametry_modelu();
    f2=(u1*T_H + u2*T_C + w1*w2 -alpha*sqrt(x1)*x2) / (C*x1^2);
end

function [y] = g(x1, x2, x2_opoz)
    y(1)=x1;
    y(2)=x2_opoz;
end  
