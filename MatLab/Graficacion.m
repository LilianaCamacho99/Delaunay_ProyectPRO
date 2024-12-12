% Leer archivo .node
fid = fopen('puntos.node', 'r');
header = textscan(fid, '%d %d %d %d', 1);
num_points = header{1};
data = textscan(fid, '%d %f %f %d', num_points);
points = [data{2} -data{3}];
fclose(fid);

% Leer archivo .ele
fid = fopen('puntos.ele', 'r');
header = textscan(fid, '%d %d %d', 1);
num_triangles = header{1};
data = textscan(fid, '%d %d %d %d', num_triangles);
triangles = double([data{2} data{3} data{4}]); % Convertir a double

% Graficar
figure;
triplot(triangles, points(:,1), points(:,2), 'b-');
hold on;
plot(points(:,1), points(:,2), 'r.', 'MarkerSize', 10);
grid on;
axis equal;
title('Triangulación de Delaunay');
xlabel('X');
ylabel('Y');