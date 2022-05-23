% Model transmitancji dla zakłóceń (stałe wartości sterowania u1 i u2)

[T_H,T_C,tau,tau_C,C0,alpha] = parametry_modelu();
Tp=1;

% dx1/dt=0.011491082919654348225776797205369*35 + 0.011491082919654348225776797205369*F_D + 0.011491082919654348225776797205369*13 - 0.0062325591795655670130200470205819*h - 0.33898110318949952112718099788515
% dx2/dt=0.010563598933309752000162527307748*35 + 0.013098862677304092480201533861608*F_D + 0.025775181397275794880396566630906*13 - 0.024929950246959761388686940810244*T + 0.0046479835306562908800715120154093*T_D - 0.0078044569166772114090629874428769*h + 0.28039692224772851125671901454039

A=[-0.0062325591795655670130200470205819, 0;
     -0.0078044569166772114090629874428769, -0.024929950246959761388686940810244];
B=[0.011491082919654348225776797205369, 0;
     0.013098862677304092480201533861608, 0.0046479835306562908800715120154093];
C = [1, 0;
     0, 1];
D = [0, 0;
     0, 0];

% Modele bez uwzględnienia opóźnień 
SS0=ss(A,B,C,D);
TF0=tf(SS0);
[num, den] = tfdata(TF0);

% Uwzględnienie opóźnień w trasmitancji
TF = tf(num, den,'IODelay',[0, 0; tau, tau])

Gz=c2d(TF,Tp,'zoh');
figure; step(Gz);
 