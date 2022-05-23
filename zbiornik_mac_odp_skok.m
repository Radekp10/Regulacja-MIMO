% Macierzowa odp. skokowa dla zbiornika

function S=zbiornik_mac_odp_skok(Tp, timefinal)
    % konstrukcja macierzowej odpowiedzi skokowej dla DMC, zakładając timefinal=Tp*D
    % Tp - okres próbkowania, D - horyzont dynamiki obiektu.
    % Model zbiornika
    s=tf('s');
    G11=0.01149/(s + 0.006233);
    G12=exp(-120*s)*0.01149/(s + 0.006233);
    G21=exp(-60*s)*(0.02578*s + 7.096e-05)/(s^2 + 0.03116*s + 0.0001554);
    G22=exp(-180*s)*(0.01056*s - 2.384e-05)/( s^2 + 0.03116*s + 0.0001554);
    G=[G11 G12;G21 G22];
       
    Gz=c2d(G,Tp,'zoh');

    % Dyskretna odpowiedź skokowa wielowymiarowa wg konwencji Matlaba:
    figure; step(Gz,timefinal); title('Odpowiedź skokowa');
    Y=step(Gz,timefinal); % macierz o wymiarze (timefinal/Tp+1) x ny x nu na odcinku
    % czasu od t=0 do t=timefinal z krokiem Tp
    % Dyskretna odpowiedź skokowa macierzowa S o wymiarze ny x nu x D, pierwsza
    % macierz dla czasu dyskr.k=1, ostatnia dla czasu dyskr. k=D=timefinal/Tp (D macierzy):
    [nt,ny,nu]=size(Y); % nt=D+1 (Y zawiera też wartość wyjść w chwili 0, macierz S nie)
    S=zeros(ny,nu,nt-1);
    for i=1:ny
        for j=1:nu
            S(i,j,:)=Y(2:nt,i,j);
        end
    end
end
