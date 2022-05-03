function result = control_error(Y, Yzad)
%% funkcja licz¹ca b³¹d œredniokwadratowy
    [sim_time, ny] = size(Y);
    result = 0;
    
    for k = 1 : sim_time
        for y_number = 1 : ny
            result = result + power(Y(k, y_number) - Yzad(k, y_number), 2);
        end
    end 
end