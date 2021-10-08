% creates a pole zero plot of the system
function fig = pole_zero_plot(sys, figname)

fig = figure('Name', figname, 'NumberTitle', 'off');
figure(fig);
fig = pzplot(sys);
xlabel('Re$(s)$')
fig.AxesGrid.XUnits = ''; 
ylabel('Im$(s)$')
fig.AxesGrid.YUnits = '';
title('')
h = findobj(gca, 'type', 'line');
set(h, 'markersize', 12) % larger markers for easier visibility
fig.AxesGrid.BackgroundAxes.XLabel.Interpreter = 'Latex'; 
fig.AxesGrid.BackgroundAxes.YLabel.Interpreter = 'Latex';

end