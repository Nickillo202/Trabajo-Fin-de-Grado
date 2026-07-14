clear all

addpath(genpath('C:/...'));               
addpath(genpath('C:/...'));      % Es donde está la función readnii (interna del grupo de investigación)
path_to_images = fullfile('C:/...');          % Ruta a las imágenes
path_to_segmentations = fullfile('C:/...');   % Ruta con las segmentaciones a corregir
path_out = fullfile('C:/...');
mkdir(path_out);
path_images_list = getFileNames(path_to_images,'*.nii.gz');
path_segmentations_list = getFileNames(path_to_segmentations,'*.nii.gz');

for i= 1:length(path_segmentations_list)
    

    full_path_to_segmentation = fullfile(path_to_segmentations,path_segmentations_list{i});
    [segmentation ,MS1,MT1]=readnii(full_path_to_segmentation);
    id = char(extractBefore(path_segmentations_list{i}, "_04"));

    image_index = find(contains(path_images_list, id));
    full_path_to_image = fullfile(path_to_images, path_images_list{image_index(1)});
    [image ,MS,MT]=readnii(full_path_to_image);

    corrected_segmentation = segmentation;         % Copia donde se guarda la corrección

    if all(segmentation(10,:,:)==3,'all') || all(segmentation(11,:,:)==3,'all') || all(segmentation(12,:,:)==3,'all') || all(segmentation(13,:,:)==3,'all') || all(segmentation(14,:,:)==3,'all') 
        segmentation_axis=1;
        segmentation_length = length(squeeze(segmentation(:,1,1)));
    elseif all(segmentation(:,10,:)==3,'all') || all(segmentation(:,11,:)==3,'all') || all(segmentation(:,12,:)==3,'all') || all(segmentation(:,13,:)==3,'all') || all(segmentation(:,14,:)==3,'all') 
        segmentation_axis=2;
        segmentation_length = length(squeeze(segmentation(1,:,1)));
    elseif all(segmentation(:,:,10)==3,'all') || all(segmentation(:,:,11)==3,'all') || all(segmentation(:,:,12)==3,'all') || all(segmentation(:,:,13)==3,'all') || all(segmentation(:,:,14)==3,'all') 
        segmentation_axis=3;
        segmentation_length = length(squeeze(segmentation(1,1,:)));
    end
        
    for slice_num= 1:5:segmentation_length
        switch segmentation_axis
            case 1
                slice_2d_segmentation = squeeze(segmentation(slice_num,:,:));
                slice_2d = squeeze(image(slice_num,:,:));
            case 2
                slice_2d_segmentation = squeeze(segmentation(:,slice_num,:));
                slice_2d = squeeze(image(:,slice_num,:));
            case 3
                slice_2d_segmentation = squeeze(segmentation(:,:,slice_num));
                slice_2d = squeeze(image(:,:,slice_num));
        end
        
        % Filtrado de imagen y segmentación
        
        filtered_slice=sobelMask(slice_2d);
        filtered_segmentation = removeCenter(slice_2d_segmentation,3);
        
        % Dejar solo píxeles diferenciales más brillantes
        threshold = quantile(filtered_slice(:),0.9);
        condition_mask = (filtered_slice > threshold);
        conditioned_filtered_slice = filtered_slice .* condition_mask;
        
        % Dejar solo píxeles más oscuros
        inverted_slice_2d=max(slice_2d(:))-slice_2d;
        threshold2 = quantile(inverted_slice_2d(:),0.9);
        condition_mask2 = (inverted_slice_2d > threshold2);
        conditioned_slice = inverted_slice_2d .* condition_mask2;
        
        % SEGMENTACIÓN CORREGIDA
        pixels_to_remove = filtered_segmentation.*~condition_mask;
        placenta_mask = slice_2d_segmentation==2;
        
        % Ahora quito la placenta de la mascara a eliminar y el contorno del fluido
        % amniótico alrededor de la placenta
        SE = strel("diamond",4);
        dilated_placenta_mask = imdilate(placenta_mask, SE);
        pixels_to_remove_no_placenta = pixels_to_remove;
        pixels_to_remove_no_placenta(dilated_placenta_mask) = 0;
        pixels_to_remove_clean = bwareaopen(pixels_to_remove_no_placenta, 5);            % Quita todos los grupos de píxeles que tengan menos de 6 píxeles unidos (unidos en diagonal cuenta)
        slice_2d_segmentation_clean = slice_2d_segmentation;
        slice_2d_segmentation_clean(pixels_to_remove_clean) = 0;

        % Guardar la segmentación corregida
        switch segmentation_axis
            case 1
                corrected_segmentation(slice_num,:,:) = slice_2d_segmentation_clean;
            case 2
                corrected_segmentation(:,slice_num,:) = slice_2d_segmentation_clean;
            case 3
                corrected_segmentation(:,:,slice_num) = slice_2d_segmentation_clean;
        end
    end
    writenii(sprintf('%s%s_04_ManualSegmentationUterusRegions02_00.nii.gz',path_out,id),corrected_segmentation,[],MS,MT);
end


% FUNCIÓN PARA APLICAR LA MÁSCARA SOBEL DE DETECCIÓN DE BORDES

function G = sobelMask(image_2d)            

    Sx = [-1, 2, -1; -1, 2, 1; -1, 2, -1];       
    Sy = Sx';

    Gx = imfilter(image_2d, Sx);
    Gy = imfilter(image_2d, Sy);
    G = sqrt(Gx.^2 + Gy.^2);
end



% FUNCIÓN PARA DEJAR SOLO EL CONTORNO DE LA SEGMENTACIÓN : Toma una
% segmentación del fluido amniótico y placenta y devuelve la segmentación 
% quitando las etiquetas de fluido amniótico que no estén aproximadamente
% a una distancia N del borde en cualquier dirección

function filtered_segmentation = removeCenter(segmentation_2d, N)
    mask = ones(1+2*N);
    filtered_segmentation_raw = imfilter(segmentation_2d, mask);
    filtered_segmentation_raw(filtered_segmentation_raw==(1+2*N)^2) = 0;
    filtered_segmentation = single(filtered_segmentation_raw & segmentation_2d);       % Las segmentaciones se pondrían guardar en uint8, pero al leer los archivos salen como single (single-precision floating-point number), así que de momento lo dejo así
    placenta_mask = segmentation_2d==2;                                                % Para restaurar la segmentación completa de la placenta, solo interesa cambiar la del fluido amniótico
    filtered_segmentation(placenta_mask) = 2;
end