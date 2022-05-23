function Sz=zbiornik_mac_odp_skok_zakl(Tp,timefinal)
    % konstrukcja macierzowej odpowiedzi skokowej zakłóceniowej dla DMC,
    % zakładając timefinal= Tp*Dz, Dz -- liczba elementów odp. skokowej
    % Model zbiornika:
    s=tf('s');
    G11=0.01149/(s + 0.006233);
    G12=0;
    G21=exp(-60*s)*(0.0131*s - 8.042e-06)/(s^2 + 0.03116*s + 0.0001554);
    G22=exp(-60*s)*(0.004648)/( s + 0.02493);
    G=[G11 G12;G21 G22];
    
    Gz=c2d(G,Tp,'zoh');
    
    % Dyskretna odpowiedź skokowa wielowymiarowa wg konwencji Matlaba:
    figure; step(Gz,timefinal); title('Odpowiedź skokowa zakłóceniowa');
    Ydstepz=step(Gz,timefinal);% macierz wymiaru (timefinal/Tp)+1 x ny x nz na odcinku
    % czasu od t=0 do t=timefinal z krokiem Tp
    % Dyskretna odpowiedź skokowa macierzowa Sz o wymiarze ny x nz x timefinal/Tp; pierwsza
    % macierz dla czasu dyskr. k=1, ostatnia dla k=Dz=timefinal/Tp (Dz macierzy):
    [nt,ny,nz]=size(Ydstepz); % nt=D+1 (Ydstepz zawiera wartość wyjść w chwili 0, Sz nie)
    Sz=zeros(ny,nz,nt-1);
    for i=1:ny
        for j=1:nz
            Sz(i,j,:)=Ydstepz(2:nt,i,j); %
        end
    end
end

