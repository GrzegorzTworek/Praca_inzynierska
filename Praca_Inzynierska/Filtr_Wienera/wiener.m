% w poni¿szej metodzie wykorzystano filtr Wienera do usuniêcia zak³óceñ.
% Metoda ta mo¿e byæ zastosowana, gdy funkcja rozmycia punktu i poziom
% zaszumienia jest znany
clc
clear all

%wczytanie obrazu
K=imread('jon-tang512.jpg');
I=rgb2gray(K);
imshow(I);
title('a) Obraz orginalny');
N=512;
%zmiana formatu obrazu na double
Obraz=double(I);
text(size(I,2),size(I,1)+15, ...
   'Photo by Jon Tang on Unsplash, link - https://unsplash.com/photos/-JODYEU6DmU', ...
   'FontSize',7,'HorizontalAlignment','right');
figure;
imhist(I);title('a) Histogram obrazu orginalnego')

%Wyznaczenie amplitudy sygna³u dla obrazu orginalnego
ObrazTF=fft2(Obraz); %przekszta³cenie obrazu do dziedzniy czêstotliwoœci(transfomacja Fouriera TF)
Amplituda_sygnalu=fftshift((abs(ObrazTF).^2)./(N*N));

%Na³o¿enie szumu na obraz o sigmie równej 40
sigma=40;
Szum=sigma*randn(size(Obraz));
figure;
imshow(uint8(Obraz+Szum));
title('a) Obraz zaszumiony');
figure;
imhist(uint8(Obraz+Szum));
title('b) Histogram obrazu zaszumionego')

%Wyznaczenie amplitudy szumu
SzumTF=fft2(Szum); % przekszta³cenie do dzidziny czêstostliwoœci(transfomacja Fouriera TF)
Amplituda_szumu=fftshift((abs(SzumTF).^2)./(N*N));
Wartosc_szumu=sum(sum(Amplituda_szumu))/(N*N)

%Transformacja Fouriera zaszumionego obrazu 
L=Obraz;
L=L+Szum;
LFT=fft2(L);
Zaszumiony_obrazTF=fftshift(LFT);

%Na³o¿enie filtru gaussowskiego na zaszumiany obraz
%Filtr gaussowski dzia³a w tym wypadku jako filtr wyg³adzaj¹cy
Filtr_wygladzajacy = fspecial('gaussian',[5 5], 1);
Wygladzenie_z_szumem = imfilter(uint8(L), Filtr_wygladzajacy, 'replicate');
figure,imshow(Wygladzenie_z_szumem);
title('a) Zaszumiony obraz po wyg³adzeniu filtrem gaussowskim');
figure;
imhist(Wygladzenie_z_szumem);
title('c) Histogram zaszumionego i wyg³adzonego obrazu')

%Amplituda sygna³u dla zaszumionego i wyg³adzonego obrazu 
Wygladzenie_z_szumemTF=fft2(Wygladzenie_z_szumem); % przekszta³cenie do dzidziny czêstostliwoœci(transfomacja Fouriera TF)
Amplituda=fftshift((abs(Wygladzenie_z_szumemTF).^2)./(N*N)); 

%Rekonstrukcja filtrem Wienera w dziedzinie czêstotliwoœci
Filtr_wypelniajacy_obraz=zeros(N); 
Filtr_wypelniajacy_obrazFT=zeros(N); 
for u=1:N 
    for v=1:N 
        %NSR (noise-to-signal response) stosunek mocy szumu do mocy sygna³u
        NSR(u,v)=((Amplituda(u,v)+Wartosc_szumu)/Amplituda(u,v)); 
        Filtr_wypelniajacy_obrazFT(u,v)=Zaszumiony_obrazTF(u,v)/NSR(u,v); 
    end
end


%odwrotna transformata Fouriera potrzebna by zwizualizowaæ obraz
OdwroconaTF_dla_obrazu=(real(ifft2(ifftshift(Filtr_wypelniajacy_obrazFT)))); 

%konwertowanie na typ uint8
figure,imshow(uint8(OdwroconaTF_dla_obrazu)); 
title('a) Odtworzony obraz filtrem Wienera')
figure;
imhist(uint8(OdwroconaTF_dla_obrazu));
title('d) Histogram obrazu odtworzonego filtrem Wienera')


%Rekonstrukcja dla NSR równego 1.5
Filtr_wypelniajacy_obraz=zeros(N);
Filtr_wypelniajacy_obrazFT=zeros(N);
for u=1:N
    for v=1:N
        %NSR (noise-to-signal response) stosunek mocy szumu do mocy sygna³u
        NSR15(u,v)=1.5*(Amplituda(u,v)/(Amplituda(u,v)+Wartosc_szumu)); 
        Filtr_wypelniajacy_obrazFT(u,v)=Zaszumiony_obrazTF(u,v)*NSR15(u,v);
    end
end


%odwrotna transformata Fouriera potrzebna by zwizualizowaæ obraz
OdwroconaTF_dla_obrazu=(real(ifft2(ifftshift(Filtr_wypelniajacy_obrazFT))));

%konwertowanie na typ uint8
figure,imshow(uint8(OdwroconaTF_dla_obrazu)); 
title('a) Odtworzony obraz filtrem Wienera dla NSR*1,5')
figure;
imhist(uint8(OdwroconaTF_dla_obrazu));
title('e) Histogram obrazu odtworzonego dla NSR*1,5')

%Rekonstrukcja dla NSR równego 0.5
Filtr_wypelniajacy_obraz=zeros(N);
Filtr_wypelniajacy_obrazFT=zeros(N);
for u=1:N
    for v=1:N
        %NSR (noise-to-signal response) stosunek mocy szumu do mocy sygna³u
        NSR05(u,v)=0.5*(Amplituda(u,v)/(Amplituda(u,v)+Wartosc_szumu)); 
        Filtr_wypelniajacy_obrazFT(u,v)=Zaszumiony_obrazTF(u,v)*NSR05(u,v);
    end
end


%odwrotna transformata Fouriera potrzebna by zwizualizowaæ obraz
OdwroconaTF_dla_obrazu=(real(ifft2(ifftshift(Filtr_wypelniajacy_obrazFT))));

%konwertowanie na typ uint8
figure,imshow(uint8(OdwroconaTF_dla_obrazu)); 
title('a) Odtworzony obraz filtrem Wienera dla NSR*0,5')
figure;
imhist(uint8(OdwroconaTF_dla_obrazu));
title('f) Histogram obrazu odtworzonego dla NSR*0,5')

