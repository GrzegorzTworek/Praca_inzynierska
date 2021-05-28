%Metoda z użyciem algorytmu Lucy-Richardsona do deblurowania obrazów. 
%Może być efektywnie stosowana, gdy znana jest funkcja rozmycia punktów FRP (operator rozmycia), ale informacja o szumie jest niewielka lub nie jest dostępna. 
%Zamazany i zaszumiony obraz zostaje przywrócony przez iteracyjny, przyspieszony, tłumiony algorytm Lucy-Richardsona. 
%Można użyć charakterystyki układu optycznego jako parametrów wejściowych, aby poprawić jakość przywracania obrazu.
 
%Funkcja dekonwolucji może obsługiwać tablice o dowolnym rozmiarz
clc
clear all
 
I = im2double(imread('jon-tang490.jpg')); 
 
figure;
imshow(I);
title('Obraz orginalny');
text(size(I,2),size(I,1)+15, ...
   'Photo by Jon Tang on Unsplash, link - https://unsplash.com/photos/-JODYEU6DmU', ...
   'FontSize',7,'HorizontalAlignment','right');
figure;
imhist(I);title('Histogram obrazu orginalnego')
 
rozmiar=[5 5]; sigma=5; 
FRP = fspecial('gaussian', rozmiar, sigma); 
Rozmycie = imfilter(I, FRP); 
figure; imshow(Rozmycie); title('Obraz rozmyty') 
figure; imhist(Rozmycie); title('Histogram rozmytego obrazu ') 
 
Srednia_szumu = 0; %wpływa na jasność
Amplituda_szumu = 0.002;
V = .002;
Rozmycie_z_szumem = imnoise(Rozmycie,'gaussian',0,V);
figure;
imshow(Rozmycie_z_szumem);
title('Rozmycie z szumem na obrazie początkowym');
figure;
imhist(Rozmycie_z_szumem);
title('Hoistogram rozmycia z szumem ');
 
tic;
wynikRLF = czestotRL(Rozmycie_z_szumem, FRP, 15); toc; 
figure; imshow(wynikRLF); title('Odtworzony obraz w dziedzinie częstotliwosci')
figure; imhist(wynikRLF); title('Histogram odtworzonego obrazu w dziedzinie częstotliwosci')
 
tic;
wynikRL = dekonwolucjaRL(Rozmycie_z_szumem, FRP, 15); toc; 
figure; imshow(wynikRL); title('Odtworzony obraz')
figure; imhist(wynikRL); title('Histogram odtworzonego obrazu')
 
%Kontrola wzmocnienia szumu poprzez tłumienie
%Najnowszy obraz, LucyRich2, jest wynikiem 15 iteracji. 
%Chociaż jest ostrzejszy niż poprzedni wynik z 5 iteracji, obraz rozwija "plamkowaty" wygląd. 
%Plamki nie odpowiadają żadnym prawdziwym strukturom (co widać poprzez porównanie z pierwotnym obrazem), ale są wynikiem zbytniego dopasowania szumu do danych.
 
 
Tlumik = im2double(4*sqrt(V));
LucyRich3 = deconvlucy(Rozmycie_z_szumem,FRP,15,Tlumik);
figure;
imshow(LucyRich3);
title('Odtworzony obraz z tłumieniem, po 15 iteracjach');
figure;
imhist(LucyRich3)
title('Histogram odtworzonego obrazu z tłumieniem, po 15 iteracjach');
 
%Zostanie teraz wykonana rekonstrukcje z funkcją Waga, która działa na samym rozmyciu obrazu. 
%Tworzymy macierz funkcji ważenia, Waga, która składa się z jedynek w środkowej części zamazanego obrazu ("dobre" piksele, znajdujące się w linii przerywanej) i zer na krawędziach ("złe" piksele - te, które nie odbierają sygnału ).
 
%Wykonano macierz Waga zgodną z romiarem obrazu orginalnego
%Algorytm waży każdą wartość piksela zgodnie z tablicą Waga podczas przywracania obrazu. 
%W poniższym przypadku używane są tylko wartości środkowych pikseli (gdzie WEIGHT = 1), podczas gdy "złe" wartości pikseli są wykluczone z optymalizacji.
 
K = im2double(Rozmycie_z_szumem)
Waga = zeros(size(I));
Waga(5:end-4,5:end-4) = 1;
CutImage = K.*Waga;
CutEdged = edgetaper(CutImage,FRP);
DAMPAR = im2double(4*sqrt(V));
LucyRich4 = deconvlucy(CutEdged,FRP,15,DAMPAR,Waga);
figure; imshow(LucyRich4)
title('Oczyszczony obraz')
figure; imhist(LucyRich4);
title('Histogram obrazu oczyszczonego z funkcja ważona')
 
 
function result = dekonwolucjaRL(obraz, FRP, iter) 
    % obrazy muszą być w formacie double by wykonać konwolucjee 
    obraz = double(obraz); 
    FRP = double(FRP); 
    obrazdys = obraz; % wprowadzenie na obraz dyskretny 
    OdwrFRP = FRP(end:-1:1,end:-1:1); % przestrzennie odwrócona PSF dla obrazu dyskretnego 
    
    for i= 1:iter 
     konwolucja  = convn(obrazdys,FRP,'same'); 
     Rozmycie2 = obraz./konwolucja; 
     blad  = convn(Rozmycie2,OdwrFRP,'same'); 
     obrazdys = obrazdys.* blad; 
    end 
    
    result = obrazdys;
end
 
%funkcja w dziedzinie czestotliwosci
 
function result = czestotRL(obraz, FRP, iter) 
obraziter = obraz; %wczytanie obrazu do iteracji 
OFP = psf2otf(FRP,size(obraz)); %OFP optyczna funkcja przenoszenia 
%z ang. OTF-optical transfer function, zmiana FRP na OTF
 
for i=1:iter 
    uklad = fft2(obraziter); %przekształcenie do dziedziny częstotliwości
    splot = OFP.*uklad; 
    transformacjaodw = ifft2(splot); %odwrotna transformacja 
    %fouriera(czestotliwosciowa)
    stosunek = obraz./transformacjaodw; 
    ukladstosunku = fft2(stosunek);  
    wynik = OFP .* ukladstosunku; 
    transfwyniku = ifft2(wynik); 
    obraziter = transfwyniku.*obraziter;  
end 
result = abs(obraziter); 
end
