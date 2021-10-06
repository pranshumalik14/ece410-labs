% automatically trim and save figure to a pdf
function savefig(fig, relfilename)

set(fig, 'Units', 'Inches');
pos = get(fig, 'Position');
set(fig, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)])
print(fig, relfilename, '-dpdf', '-r0')