function y = DownSample(x, n, phase)
% DOWNSAMPLE Réduction du taux d'échantillonnage d'un facteur entier
%   y = downsample(x, n) réduit le taux d'échantillonnage de x en gardant
%   un échantillon sur n. Le premier échantillon est toujours conservé.
%
%   y = downsample(x, n, phase) spécifie un décalage de phase (offset)
%   pour commencer l'échantillonnage. phase doit être un entier entre 0 et n-1.
%
%   Exemples:
%       x = [1 2 3 4 5 6 7 8 9 10];
%       y = downsample(x, 3)      % Résultat: [1 4 7 10]
%       y = downsample(x, 3, 1)   % Résultat: [2 5 8]
%
%   Voir aussi upsample, decimate, resample

% Validation des arguments d'entrée
if nargin < 2 || nargin > 3
    error('Nombre incorrect d''arguments. Utilisation: downsample(x, n) ou downsample(x, n, phase)');
end

if nargin == 2
    phase = 0; % Valeur par défaut
end

% Vérification des paramètres
if ~isscalar(n) || n ~= floor(n) || n <= 0
    error('Le facteur n doit être un entier positif');
end

if ~isscalar(phase) || phase ~= floor(phase) || phase < 0 || phase >= n
    error('Le paramètre phase doit être un entier entre 0 et n-1');
end

% Conversion en colonne si x est un vecteur ligne
if isrow(x)
    x = x.';
    row_output = true;
else
    row_output = false;
end

% Sélection des échantillons
y = x(phase+1:n:end, :);

% Conversion en ligne si l'entrée était un vecteur ligne
if row_output && size(y, 1) == 1
    y = y.';
elseif row_output
    y = y.';
end
end