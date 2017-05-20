Scale = [0.424,0.847,0.558,0.000,0.597,2.007,0.000,0.913,2.353,0.231,0.641,0.000,0.000,0.481,0.341];
Shape = [0.215,0.000,0.157,0.344,0.391,0.000,0.500,0.750,0.750,0.318,0.122,0.393,0.421,0.000,0.157];
Name = {' Equirectangular',sprintf('Mercator\n'),' Gall stereographic',' Hobo-Dyer',' Polar', sprintf('Stereographic\n'),' Azimuthal equal-area',' Orthographic','Gnomonic ',' Winkel tripel',' Van der Grinten',' Mollweide',' Hammer',sprintf('Pierce quincuncial\n'),'Tetragraph '};
Align = {'left','left','left','left','left','left','left','left','right','left','left','left','left','right','right'};

plot(Scale, Shape, '.', 'MarkerSize',20);
xlabel('Size distortion');
ylabel('Shape distortion');
title('Map projections, arranged by distortion');
for i = 1:length(Scale)
    text(Scale(i), Shape(i), Name{i}, 'HorizontalAlignment', Align{i});
end