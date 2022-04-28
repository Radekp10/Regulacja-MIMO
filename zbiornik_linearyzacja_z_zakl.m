% Automatyczna linearyzacja przy pomocy pakietu symbolicznego Matlaba
syms F_C F_H h T F_D T_D

F_C0 = 35;
F_H0 = 13;
F_D0 = 11;
T_D0 = 31;
h0 = 54.39;
T0 = 34.05;

C = 0.8;
a = 8;
T_C = 25;
T_H = 61;
tau_C = 120;
tau = 60;

f1(F_C, F_H, h, T, F_D, T_D) = 1 / (2 * C * h) * (F_H + F_C + F_D - a * sqrt(h));
f2(F_C, F_H, h, T, F_D, T_D) = 1 / (C * h^2) * (F_H * T_H + F_C * T_C + F_D * T_D - a * sqrt(h) * T);

% Prawdopodobnie idzie ³adniej, ale tak te¿ dzia³a (nie jest odporny na
% zmianê punktu pracy !!!!!!)
feval(symengine, 'mtaylor', f1, '[h = 54.39, T = 34.05, F_C=35, F_H=13, F_D=11, T_D=31]', 2)
feval(symengine, 'mtaylor', f2, '[h = 54.39, T = 34.05, F_C=35, F_H=13, F_D=11, T_D=31]', 2)
