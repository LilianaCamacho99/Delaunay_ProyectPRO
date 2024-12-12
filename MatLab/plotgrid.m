% Nombre base del archivo
basename = 'puntos';

% Cargar archivo .node
nodefile = [basename '.node'];
fid = fopen(nodefile);                        % Cargar archivo de conectividad basado en vértices de TRIANGLE
[nnode] = fscanf(fid, '%i', [1 4]);           % Obtener número de nodos
ncol = 3 + nnode(3) + nnode(4);               % Especificar número de columnas en el archivo .node
data = fscanf(fid, '%f', [ncol nnode(1)])';   % Obtener datos
x = data(:, 2); 
y = -data(:, 3);                              % Reflejar las coordenadas del eje y
fclose(fid);

% Cargar archivo .ele
elefile = [basename '.ele'];
fid = fopen(elefile);                         % Cargar archivo de conectividad basado en elementos de TRIANGLE
[nelem] = fscanf(fid, '%i', [1 3]);           % Obtener número de triángulos
ncol = 4 + nelem(3);                          % Especificar número de columnas en el archivo .ele
tri = fscanf(fid, '%i', [ncol nelem(1)])';    % Obtener tabla de conectividad
fclose(fid);

% Crear la malla usando trimesh
trimesh(tri(:, 2:4), x, y, zeros(size(x)), ...   
                  'EdgeColor', 'k', ...
                  'FaceColor', 'none', ...
                  'LineWidth', 0.5);

% Configurar vista 2D
view(2);                                       % Vista en 2D
axis equal;                                    % Ejes iguales
