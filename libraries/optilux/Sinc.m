function y = Sinc(x)
% SINC Fonction sinus cardinal
%   y = sinc(x) calcule le sinus cardinal de x, défini comme :
%   sinc(x) = sin(πx)/(πx) pour x ≠ 0
%   sinc(0) = 1
%
%   La fonction fonctionne pour les scalaires, vecteurs, matrices et tableaux N-D
%
%   Exemples:
%       sinc(0)   % Retourne 1
%       sinc(1)   % Retourne 0
%       sinc([-2 -1 0 1 2]) % Retourne [0 0 1 0 0]
%
%   Voir aussi sin, diric, gauspuls

% Calcul du sinus cardinal
y = ones(size(x), 'like', x); % Initialisation avec des 1 (pour gérer x=0)
i = (x ~= 0);                % Indices où x n'est pas nul

% Calcul pour les éléments non nuls
if any(i(:))
    x_nonzero = x(i);
    y(i) = sin(pi*x_nonzero) ./ (pi*x_nonzero);
end

% Gestion des nombres complexes
if ~isreal(x)
    y = complex(y);
end
end