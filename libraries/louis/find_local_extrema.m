function [minima, maxima] = find_local_extrema(Lkm, LX, threshold)
% Trouve les abscisses des extrema locaux (minimums et maximums) d'une courbe LX = f(Lkm).
% Entrées :
%   Lkm : Vecteur des abscisses (1, 2, ..., 200)
%   LX : Vecteur des ordonnées (valeurs de LX pour chaque Lkm)
%   threshold : Seuil pour filtrer les petites variations (par ex., 1e-4)
% Sorties :
%   minima : Vecteur des abscisses des minimums locaux
%   maxima : Vecteur des abscisses des maximums locaux

% Vérification des entrées
if length(Lkm) ~= length(LX)
    error('Les vecteurs Lkm et LX doivent avoir la même longueur.');
end
if length(Lkm) < 3
    error('Il faut au moins 3 points pour détecter des extrema.');
end

% Initialisation des sorties
minima = [];
maxima = [];

% Calcul de la dérivée numérique
dLX = diff(LX) ./ diff(Lkm);

% Boucle pour trouver les changements de signe dans la dérivée
for i = 1:(length(dLX)-1)
    % Produit des dérivées consécutives
    product = dLX(i) * dLX(i+1);
    
    % Vérifier si la dérivée change de signe (produit négatif)
    if product < 0
        % Estimer l'abscisse de l'extremum comme le point médian
        x_extremum = Lkm(i+1);
        
        % Vérifier si la variation est significative (supérieure au seuil)
        if abs(dLX(i)) > threshold && abs(dLX(i+1)) > threshold
            % Déterminer si c'est un minimum ou un maximum
            % Si dLX(i) > 0 et dLX(i+1) < 0, c'est un maximum
            if dLX(i) > 0
                maxima = [maxima; x_extremum];
            % Si dLX(i) < 0 et dLX(i+1) > 0, c'est un minimum
            elseif dLX(i) < 0
                minima = [minima; x_extremum];
            end
        end
    end
end

% Affidisplay('Abscisses des minimums locaux :', minima);
% display('Abscisses des maximums locaux :', maxima);
end