function plot_RAP_txmap(figure,s,headline,color_range,color_data,fpos)
% modify a plot for tx map
set(figure,'Units','normalized','Position',[fpos 0.4 0.5])
set(s,'Marker','o','MarkerFaceColor','flat','SizeData',10,'CData',color_data*1000)

xlim([-158.3 -157.7])
ylim([22.45 23])
grid on
xlabel('Longitude')
ylabel('Latitude')



cbar = colorbar;
caxis([color_range]);
colormap jet;
cbar.Label.String = 'Travel Time Perturbation (ms)';

set(gca,'fontsize',16)
title(headline,'fontsize',18)






