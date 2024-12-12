% Cargar y mostrar la imagen
imagen = imread('acapulco.jpg');
imshow(imagen);
title('Haz clic para seleccionar puntos de la Región 1. Presiona ESPACIO para Región 2, ENTER para finalizar.');

% Inicializar arrays para ambas regiones
x1 = []; y1 = [];
x2 = []; y2 = [];
region1_complete = false;

% Recolectar puntos para ambas regiones'
while true
    [x, y, button] = ginput(1);
    
    if isempty(button) % Si se presiona Enter
        break;
    elseif button == 32 % Si se presiona Espacio
        region1_complete = true;
        title('Selecciona puntos de la Región 2. Presiona ENTER para finalizar.');
        continue;
    end
    
    if ~region1_complete
        x1 = [x1; x];
        y1 = [y1; y];
        hold on;
        plot(x, y, 'r.', 'MarkerSize', 10);
        if length(x1) > 1
            plot([x1(end-1), x1(end)], [y1(end-1), y1(end)], 'r-');
        end
    else
        x2 = [x2; x];
        y2 = [y2; y];
        hold on;
        plot(x, y, 'b.', 'MarkerSize', 10);
        if length(x2) > 1
            plot([x2(end-1), x2(end)], [y2(end-1), y2(end)], 'b-');
        end
    end
end

% Cerrar los polígonos visualmente
if ~isempty(x1)
    plot([x1(end), x1(1)], [y1(end), y1(1)], 'r-');
end
if ~isempty(x2)
    plot([x2(end), x2(1)], [y2(end), y2(1)], 'b-');
end

% Crear nodos para ambas regiones
nodes1 = [(1:length(x1))', x1, y1, ones(length(x1), 1)];
nodes2 = [(length(x1)+1:length(x1)+length(x2))', x2, y2, 2*ones(length(x2), 1)];
nodes = [nodes1; nodes2];

% Generar segmentos para ambas regiones
segments1 = [(1:size(nodes1, 1))', [2:size(nodes1, 1), 1]', (1:size(nodes1, 1))', ones(size(nodes1, 1), 1)];
if ~isempty(nodes2)
    n1 = size(nodes1, 1);
    segments2 = [(n1+1:size(nodes, 1))', [n1+2:size(nodes, 1), n1+1]', (n1+1:size(nodes, 1))', 2*ones(size(nodes2, 1), 1)];
    segments = [segments1; segments2];
else
    segments = segments1;
end

% Abrir el archivo para escritura
fileID = fopen('puntos.poly', 'w');

% Escribir nodos
fprintf(fileID, '%d 2 0 1\n', size(nodes, 1));
for i = 1:size(nodes, 1)
    fprintf(fileID, '%d %.4f %.4f %d\n', nodes(i, 1), nodes(i, 2), nodes(i, 3), nodes(i, 4));
end

% Escribir segmentos
fprintf(fileID, '%d 1\n', size(segments, 1));
for i = 1:size(segments, 1)
    fprintf(fileID, '%d %d %d %d\n', segments(i, :));
end

% Escribir huecos
fprintf(fileID, '1\n');
fprintf(fileID, '1 %.2f %.2f\n', mean(x2), mean(y2));
fprintf(fileID, '0\n');

fclose(fileID);

disp('Archivo .poly generado con éxito.');