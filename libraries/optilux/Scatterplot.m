function varargout = Scatterplot(x, n, offset, plotstring, scatfig)
% SCATTERPLOT Version finale corrigée - Gère tous les cas d'utilisation

% Vérification des arguments d'entrée
if nargin < 1
    error('Pas assez d''arguments d''entrée');
end

% Valeurs par défaut
if nargin < 2 || isempty(n)
    n = 1;
end

if nargin < 3 || isempty(offset)
    offset = 0;
end

% Gestion de la figure avant de vérifier hold
if nargin < 5 || isempty(scatfig)
    scatfig = figure;
    hold_state = false;
else
    figure(scatfig);
    hold_state = ishold;
end

% Gestion automatique des couleurs et styles
if nargin < 4 || isempty(plotstring)
    if hold_state
        % Mode hold on - utilisation du ColorOrder
        ax = gca;
        color_order = ax.ColorOrder;
        next_color = color_order(mod(ax.ColorOrderIndex-1, size(color_order,1))+1,:);
        marker_color = next_color;
        plotstring = '.';
    else
        % Mode normal - bleu par défaut
        marker_color = [0 0.4470 0.7410]; % Bleu MATLAB par défaut
        plotstring = '.';
    end
else
    % Si plotstring est spécifié
    if ~isempty(plotstring)
        marker_color = plotstring(1);
        if length(plotstring) > 1
            plotstring = plotstring(2:end);
        else
            plotstring = '.';
        end
    else
        marker_color = 'b';
        plotstring = '.';
    end
end

% Validation des paramètres
% validateattributes(x, {'numeric'}, {'vector'}, 'scatterplot', 'x', 1);
% validateattributes(n, {'numeric'}, {'positive', 'integer', 'scalar'}, 'scatterplot', 'n', 2);
% validateattributes(offset, {'numeric'}, {'nonnegative', 'integer', 'scalar', '<', n}, 'scatterplot', 'offset', 3);

% Sous-échantillonnage
x_sampled = x(offset+1:n:end);
real_part = real(x_sampled);
imag_part = imag(x_sampled);

% Tracé principal
if ischar(marker_color)
    h = plot(real_part, imag_part, [marker_color plotstring], 'MarkerSize', 8);
else
    h = plot(real_part, imag_part, plotstring, 'MarkerSize', 8, 'Color', marker_color);
end
hold on;

% Configuration des axes
grid on;
axis equal;
xlabel('Partie Réelle');
ylabel('Partie Imaginaire');
title('Diagramme de Dispersion');

% Lignes de référence
plot([0 0], ylim, 'k-', 'LineWidth', 0.5);
plot(xlim, [0 0], 'k-', 'LineWidth', 0.5);

% Ajustement des limites
current_xlim = xlim;
current_ylim = ylim;
new_xlim = [min(min(real_part), current_xlim(1)), max(max(real_part), current_xlim(2))];
new_ylim = [min(min(imag_part), current_ylim(1)), max(max(imag_part), current_ylim(2))];
max_limit = max(max(abs(new_xlim)), max(abs(new_ylim))) * 1.1;
if ~isempty(max_limit) && isfinite(max_limit)
    xlim([-max_limit max_limit]);
    ylim([-max_limit max_limit]);
end

if ~hold_state
    hold off;
end

% Retourne le handle si demandé
if nargout > 0
    varargout{1} = h;
end
end