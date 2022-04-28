% Analiza RGA obiektu

K=[0.01149/0.006233 0.01149/0.006233; 0.00007096/0.0001554 -0.00002384/0.0001554];
K1=(K^-1)';
Lambda=K.*K1

% Sugerowane połączenia: 1-2 (y1 z u2) i 2-1 (y2 z u1)