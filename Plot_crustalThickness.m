% Parameters
px = 10; % pixels per degree
filename = ['data/thick-2-38.0-8-3300.txt'];
Mx = txt2mx(filename);
% Plot
figure;
imagesc([-180 180], -[-90 90], Mx);
set(gca, 'YDir', 'normal');  % so latitude increases upward
demcmap(Mx, 128); % Using a better colormap for topography visualization
colorbar();

%title('Crustal Thickness');
xlabel('Longitude (deg)');
ylabel('Latitude (deg)');
