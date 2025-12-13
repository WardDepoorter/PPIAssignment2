% Parameters
px = 4; % pixels per degree
filename = ['Topography/LDEM_4.IMG'];
Mx = img2mx(filename, px);
Mx = [Mx(:, 181*px:end), Mx(:, 1:180*px)]; % Shift the center of the map 180 degrees

% Plot
figure;
imagesc([-180 180], -[-90 90], Mx);
set(gca, 'YDir', 'normal');  % so latitude increases upward
demcmap(Mx); % Using a better colormap for topography visualization
colorbar();

title('Lunar topography in km from reference speriod');
xlabel('Longitude (deg)');
ylabel('Latitude (deg)');
