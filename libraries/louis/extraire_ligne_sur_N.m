function extraire_ligne_sur_N(fichier_entree, fichier_sortie, N)
    % Ouvrir le fichier d'entrée en lecture
    fid_in = fopen(fichier_entree, 'r');
    if fid_in == -1
        error('Impossible d''ouvrir le fichier d''entrée');
    end
    
    % Ouvrir le fichier de sortie en écriture
    fid_out = fopen(fichier_sortie, 'w');
    if fid_out == -1
        fclose(fid_in);
        error('Impossible d''ouvrir le fichier de sortie');
    end
    
    ligne_num = 0;
    while ~feof(fid_in)
        ligne = fgetl(fid_in);
        ligne_num = ligne_num + 1;
        
        % Écrire une ligne sur N
        if mod(ligne_num - 1, N) == 0
            fprintf(fid_out, '%s\n', ligne);
        end
    end
    
    % Fermer les fichiers
    fclose(fid_in);
    fclose(fid_out);
    
    fprintf('Extraction terminée. Fichier créé : %s\n', fichier_sortie);
end