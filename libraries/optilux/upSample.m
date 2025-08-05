function y = upsample(x, L, phase)
% UPSAMPLE Augmente le taux d'échantillonnage d'un facteur entier
%   y = upsample(x, L) augmente le taux d'échantillonnage de x par un facteur L
%   en insérant L-1 zéros entre chaque échantillon. Le premier échantillon de x
%   correspond au premier échantillon de y.
%
%   y = upsample(x, L, phase) spécifie un décalage de phase (offset) pour
%   l'insertion des zéros. phase doit être un entier entre 0 et L-1.
%
%   Exemple:
%       x = [1 2 3 4];
%       y = upsample(x, 3);      % Résultat: [1 0 0 2 0 0 3 0 0 4 0 0]
%       y = upsample(x, 3, 1);   % Résultat: [0 1 0 0 2 0 0 3 0 0 4 0]
%
%   Voir aussi downsample, interp, resample

% Vérification des arguments d'entrée
if nargin < 2 || nargin > 3
    error('Nombre incorrect d''arguments d''entrée. Utilisation: upsample(x, L) ou upsample(x, L, phase)');
end

if nargin == 2
    phase = 0; % Valeur par défaut si phase n'est pas spécifiée
end

% Validation des paramètres
if ~isscalar(L) || L ~= floor(L) || L <= 0
    error('Le facteur d''upsampling L doit être un entier positif');
end

if ~isscalar(phase) || phase ~= floor(phase) || phase < 0 || phase >= L
    error('Le paramètre phase doit être un entier entre 0 et L-1');
end

% Conversion en colonne si x est un vecteur ligne
if isrow(x)
    x = x.';
    row_output = true;
else
    row_output = false;
end

% Taille de l'entrée
[m, n] = size(x);

% Initialisation de la sortie avec des zéros
y = zeros(m*L + max(0, phase), n, class(x));

% Insertion des échantillons avec le décalage approprié
y(phase+1:L:end, :) = x;

% Conversion en ligne si l'entrée était un vecteur ligne
if row_output && size(y, 1) == 1
    y = y.';
elseif row_output
    y = y.';
end
end