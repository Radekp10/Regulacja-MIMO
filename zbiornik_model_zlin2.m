% Zlinearyzowany model zbiornika (tym razem również względem zakłóceń)

function [f1Lh, f2Lh, gLh] = zbiornik_model_zlin2
    f1Lh=@f1L;
    f2Lh=@f2L;
    gLh=@gL;
end


% Zbiornik z mieszaniem  - model w przestrzeni stanu: dynamiczny, ciągły i zlinearyzowany
function [f1] = f1L(u1, u2, w1, w2, x1, x2)
    % Równanie zlinearyzowane
    f1=0.011491082919654348225776797205369*u2 + 0.011491082919654348225776797205369*w1 + 0.011491082919654348225776797205369*u1 - 0.0062325591795655670130200470205819*x1 - 0.33898110318949952112718099788515;
end

function [f2] = f2L(u1, u2, w1, w2, x1, x2)
    % Równanie zlinearyzowane
    f2=0.010563598933309752000162527307748*u2 + 0.013098862677304092480201533861608*w1 + 0.025775181397275794880396566630906*u1 - 0.024929950246959761388686940810244*x2 + 0.0046479835306562908800715120154093*w2 - 0.0078044569166772114090629874428769*x1 + 0.28039692224772851125671901454039;
end

function [y] = gL(x1, x2, x2_opoz)
    y(1)=x1;
    y(2)=x2_opoz;
end   
