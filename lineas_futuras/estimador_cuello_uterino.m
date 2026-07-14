clear all; close all;

addpath(genpath('C:/Users/nicki/Documents/TFG/NNU/myCode'));
addpath(genpath('C:/Users/nicki/Documents/TFG/Nicolas/CODE/Matlab'));
addpath(genpath('C:/Users/nicki/Documents/TFG/Nicolas/CERVIX_INTERNAL_OS_LOCATIONS'));
path_to_predictions = fullfile('C:/Users/nicki/Documents/TFG/Nicolas/EXP_2_POSTPROCESSED_A_B_PREDICTIONS_TERCER_ENTRENAMIENTO/');

% JSON donde se guardan las localizaciones identificadas manualmente
fid = fopen('Discrete_locations.json');
raw = fread(fid, inf);
str = char(raw');
fclose(fid);

discrete_locations = jsondecode(str);
case_fields        = fieldnames(discrete_locations);
num_cases          = length(case_fields);

results_CaseID   = strings(num_cases, 1);
results_DistMaxD = zeros(num_cases, 1);


for c = 1:num_cases

    matlab_case_id = case_fields{c};
    case_id        = regexprep(matlab_case_id, '^x|X', '');

    target_loc = discrete_locations.(matlab_case_id);
    target_loc = target_loc(:)';
    target_loc = [target_loc(2), target_loc(1), target_loc(3)];

    filename             = sprintf('test_%s.nii.gz', case_id);
    case_prediction_path = fullfile(path_to_predictions, filename);

    if ~isfile(case_prediction_path)
        fprintf('Warning: File not found for case %s. Skipping...\n', case_id);
        continue;
    end

    [case_prediction_file, ~, ~] = readnii(case_prediction_path);


    mask_AF       = (case_prediction_file == 1);
    mask_placenta = (case_prediction_file == 2);
    mask_body     = (case_prediction_file == 3);

    [y_all, x_all, z_all] = ind2sub(size(case_prediction_file), find(case_prediction_file > 0));
    coords_full = [x_all, y_all, z_all];
    centroid    = mean(coords_full, 1);

    % solo se considera la parte de la segmentación inferior al centroide
    coords_inf  = coords_full(z_all <= centroid(3), :);


    % Distancia entre centroide y punto de la segmentación más lejano

    dists   = sqrt(sum((coords_inf - centroid).^2, 2));
    [~, idx] = max(dists);

    P1       = centroid;
    P2       = coords_inf(idx, :);
    MaxD_inf = P2;


    % Distancia a la localización marcada
    dist_maxd        = norm(target_loc - MaxD_inf);
    results_CaseID(c)   = case_id;
    results_DistMaxD(c) = dist_maxd;


    % PLOT
    figure('Name', sprintf('3D Render: Case %s', case_id), 'Color', 'w');
    hold on;

    p1 = patch(isosurface(smooth3(mask_AF,       'box', 5), 0.5));
    set(p1, 'FaceColor', 'cyan',  'EdgeColor', 'none', 'FaceAlpha', 0.1);
    p2 = patch(isosurface(smooth3(mask_placenta, 'box', 5), 0.5));
    set(p2, 'FaceColor', 'red',   'EdgeColor', 'none', 'FaceAlpha', 0.5);
    p3 = patch(isosurface(smooth3(mask_body,     'box', 5), 0.5));
    set(p3, 'FaceColor', 'green', 'EdgeColor', 'none', 'FaceAlpha', 0.8);

    p_maxd = plot3([P1(1) P2(1)], [P1(2) P2(2)], [P1(3) P2(3)], ...
                   'k--', 'LineWidth', 3);

    xl = xlim; yl = ylim;
    fill3([xl(1) xl(2) xl(2) xl(1)], [yl(1) yl(1) yl(2) yl(2)], ...
          repmat(centroid(3), 1, 4), ...
          'yellow', 'FaceAlpha', 0.10, 'EdgeColor', 'yellow', 'LineStyle', '--');

    plot3(centroid(1), centroid(2), centroid(3), 'ys', 'MarkerSize', 12, 'MarkerFaceColor', 'yellow', 'LineWidth', 2);
    plot3(MaxD_inf(1), MaxD_inf(2), MaxD_inf(3), 'k*', 'MarkerSize', 12, 'LineWidth', 2);
    p_loc = plot3(target_loc(1), target_loc(2), target_loc(3), 'bp', 'MarkerSize', 15, 'MarkerFaceColor', 'b');

    view(3); axis tight equal; grid on;
    camlight('headlight'); lighting gouraud;
    xlabel('X'); ylabel('Y'); zlabel('Z (Slices)');
    title(sprintf('Case: %s | Dist MaxD: %.2f', case_id, dist_maxd), 'Interpreter', 'none');
    legend([p1, p2, p3, p_maxd, p_loc], {'AF', 'Placenta', 'Body', 'Max Distance Axis', 'JSON Discrete Loc'}, 'Location', 'best');
    hold off;
end


% Tabla con resultados de los 6 casos
valid_idx        = results_CaseID ~= "";
results_CaseID   = results_CaseID(valid_idx);
results_DistMaxD = results_DistMaxD(valid_idx);

ResultsTable = table(results_CaseID, results_DistMaxD, 'VariableNames', {'CaseID', 'Dist_to_MaxDiam_Inf'});

disp(' ');
disp('=================================================================');
disp('                 TABLA RESULTADOS (Vóxeles)                      ');
disp('=================================================================');
disp(ResultsTable);
fprintf('\nAverage Distance to Max Distance Inferior Pole: %.2f voxels\n', mean(results_DistMaxD));