% Parameters
filename = ['data/thick-2-38.0-8-3300.txt']; %from K25 datasets
Mx = txt2mx(filename);
data.vec.tot = sqrt(data.vec.X.^2 + data.vec.Y.^2+data.vec.Z.^2);

% Plot
figure;
imagesc([-180 180], [-90 90], data.vec.tot);
set(gca, 'YDir', 'normal');  % so latitude increases upward
demcmap(data.vec.tot, 128); % Using a better colormap for topography visualization
colorbar();

%title('Crustal Thickness');
xlabel('Longitude (deg)');
ylabel('Latitude (deg)');
data.vec.tot = sqrt(data.vec.X.^2 + data.vec.Y.^2+data.vec.Z.^2);
