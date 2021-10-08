function pole_zero_plot(sys)

figure;
f2 = pzplot(sys);

xlabel('Re$(s)$')
f2.AxesGrid.XUnits = ''; 
ylabel('Im$(s)$')
f2.AxesGrid.YUnits = '';
title('')

h = findobj(gca, 'type', 'line');
set(h, 'markersize', 15)
f2.AxesGrid.BackgroundAxes.XLabel.Interpreter = 'Latex'; 
f2.AxesGrid.BackgroundAxes.YLabel.Interpreter = 'Latex';

grid on;
axis equal;

end