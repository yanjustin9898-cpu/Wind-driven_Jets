figure('Position', [100 100 2200 2000]); % Recommended wider aspect ratio for 1:2:1 columns
indsp = 1:120;
adv_regionp = squeeze(nanmean(advectionp_highpassed(:,indsp,:),2));
sst_regionp = squeeze(nanmean(sstmean2p(:,indsp,:),2));

indsa = 201:321;
adv_regiona = squeeze(nanmean(advectiona_highpassed(:,indsa,:),2));
sst_regiona = squeeze(nanmean(sstmean2a(:,indsa,:),2));

% -------------------------------------------------------------------------
% Grid Layout Geometry (2 rows x 3 columns, Width Ratio 1:2:1)
% -------------------------------------------------------------------------
marg_h = [0.08, 0.08]; % [bottom, top]
marg_w = [0.04, 0.04]; % [left, right]
gap_h  = 0.07;         % Vertical gap
gap_w  = [0.05, 0.05]; % Horizontal gaps [Col 1-2, Col 2-3]

axh = (1 - sum(marg_h) - gap_h) / 2;

total_avail_w = 1 - sum(marg_w) - sum(gap_w);
axw1 = total_avail_w * (1.3 / 4);
axw2 = total_avail_w * (1.7 / 4);
axw3 = total_avail_w * (1 / 4);

px = [ ...
    marg_w(1), ...
    marg_w(1) + axw1 + gap_w(1), ...
    marg_w(1) + axw1 + gap_w(1) + axw2 + gap_w(2) ...
];

py = [ ...
    1 - marg_h(2) - axh, ...
    marg_h(1) ...
];

ss = zeros(1, 6);

% -------------------------------------------------------------------------
% Subplot 1 (Row 1, Col 1): Pacific Heatfluxes & SST Contour
% -------------------------------------------------------------------------
ss(1) = subplot('Position', [px(1) py(1) axw1 axh]);
contourf(longitude_p, latitude_p, advectionmeanp_highpassed, -5e-7:5e-9:5e-7, 'LineStyle', 'None'); hold on;
contour(longitude_p, latitude_p, sstp_mean_highpassed, 0:100:100, 'LineColor', 'black', 'LineWidth', 1.2);
contour(longitude_p, latitude_p, sstp_mean_highpassed, [-1e-2,0], 'LineColor', 'black', 'LineStyle', '--');
contour(longitude_p, latitude_p, sstp_mean_highpassed, [0 1e-2], 'LineColor', 'black');
cm = redblue(101); colormap(gca, cm); worldmap3(2);
cb = colorbar();  cb.Ticks = -8e-8:4e-8:8e-8; set(cb, 'AxisLocation', 'out', 'FontSize', 14);
xlim([210 290]); ylim([-30 30]); caxis([-1e-7,1e-7]);
xticks([210 240 270]); xticklabels({'150°W','120°W','90°W'});
yticks([flip(-1*[0 10 20 30]) [10 20 30]]); yticklabels({'30°S','20°S','10°S','0°','10°N','20°N','30°N'});
ax = gca; set(ax, 'FontSize', 15, 'Color', 0.6*[1 1 1], 'TickDir', 'out', 'DataAspectRatioMode', 'auto', 'PlotBoxAspectRatioMode', 'auto');
title('a)', 'Units', 'normalized', 'Position', [0, 1.02], 'HorizontalAlignment', 'left', 'FontSize', 15);

% -------------------------------------------------------------------------
% Subplot 2 (Row 1, Col 2 - Wide): Atlantic Heatfluxes & SST Contour
% -------------------------------------------------------------------------
ss(2) = subplot('Position', [px(2) py(1) axw2 axh]);

contourf(1:1461,latitude_p,adv_regionp,-1e-7:1e-9:1e-7,'edgecolor','none');
hold on;
contour(1:1461,latitude_p,sst_regionp, [-1e-2, 0], 'LineColor', 'black', 'LineStyle', '--');
contour(1:1461,latitude_p,sst_regionp, 0:100:100, 'LineColor', 'black', 'LineWidth', 1.2);
contour(1:1461,latitude_p,sst_regionp, [0 1e-2], 'LineColor', 'black');
cm = redblue(101);
colormap(cm);
cb = colorbar();
cb.Ticks = -4e-8:2e-8:4e-8; set(cb, 'AxisLocation', 'out', 'FontSize', 14);
ylim([-30 30]);
caxis([-5e-8,5e-8])
yticks([flip(-1*[0 10 20 30]) [10 20 30]]);
yticklabels({'30°S','20°S','10°S','0°','10°N','20°N','30°N'});
xticks([44,44+365,44+365*2,45+365*3]);
xticklabels({'2006','2007','2008','2009'});
ax = gca; set(ax, 'FontSize', 15);
title('b)','Units','normalized','Position',[0, 1],'HorizontalAlignment','left','FontSize',15);

% -------------------------------------------------------------------------
% Subplot 3 (Row 1, Col 3): Pacific Virtual Pot Temp Grad Contour
% -------------------------------------------------------------------------
m = size(adv_regionp, 1);
n = size(adv_regionp,2);
r = zeros(m, 1);
p = zeros(m, 1);

for i = 1:m
    [R, P] = corrcoef(sst_regionp(i,:), adv_regionp(i,:));
    r(i) = R(1, 2);
    p(i) = P(1, 2);
end

% --- 底部 x 轴（相关系数）---
ax1 = axes('Position', [px(3) py(1) axw3 axh], 'XAxisLocation', 'bottom', 'YAxisLocation', 'left');
plot(ax1, movmean(r,12), latitude_p, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 4);
%xlabel(ax1, 'Correlation coefficient (r)', 'FontSize', 15);
ylim([-30 30]);
yticks([flip(-1*[0 10 20 30]) [10 20 30]]);
yticklabels({'30°S','20°S','10°S','0°','10°N','20°N','30°N'});
ax1.XColor = 'b';
ax1.YColor = 'k';
grid(ax1, 'on');
xlim(ax1, [-1.1, 1.1]);
set(ax1, 'FontSize', 15);

% --- 顶部 x 轴（p 值）---
ax2 = axes('Position', [px(3) py(1) axw3 axh], 'XAxisLocation', 'top', 'YAxisLocation', 'left', ...
           'Color', 'none', 'Box', 'off');
ax2.YAxis.Visible = 'off';
ax2.YLim = ax1.YLim;
hold(ax2, 'on');
plot(ax2, p, latitude_p, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 4);
xlabel(ax2, 'p-value', 'FontSize', 15);
ax2.XColor = 'r';
xlim(ax2, [-0.15, 2.05]);
set(ax2, 'FontSize', 15);

xline(ax2, 0.05, 'k--', 'LineWidth', 1.5);
text(ax2, 0.05, double(max(latitude_p)*0.9), ' p = 0.05', 'Color', 'k', ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 15);

h1 = findobj(ax1, 'Type', 'line');
h2 = findobj(ax2, 'Type', 'line');
lgd = legend([h1(1), h2(1)], {'r', 'p-value'}, ...
    'Position', [px(3) + 0.035, py(1) + 0.2, 0.04, 0.03], 'FontSize', 13);


title('c)', 'Units', 'normalized', 'Position', [0, 1.02], 'HorizontalAlignment', 'left', 'FontSize', 14);

% -------------------------------------------------------------------------
% Subplot 4 (Row 2, Col 1): Pacific Heatfluxes Binned Plot
% -------------------------------------------------------------------------
ss(4) = subplot('Position', [px(1) py(2) axw1 axh]);
contourf(longitude_a, latitude_a, advectionmeana_highpassed, -5e-7:5e-9:5e-7, 'LineStyle', 'None'); hold on;
contour(longitude_a, latitude_a, ssta_mean_highpassed, 0:100:100, 'LineColor', 'black', 'LineWidth', 1.2);
contour(longitude_a, latitude_a, ssta_mean_highpassed, [-1e-2,0], 'LineColor', 'black', 'LineStyle', '--');
contour(longitude_a, latitude_a, ssta_mean_highpassed, [0 1e-2], 'LineColor', 'black');
cm = redblue(101); colormap(gca, cm); worldmap3Atl(2);
cb = colorbar(); cb.Ticks = -8e-8:4e-8:8e-8; set(cb, 'AxisLocation', 'out', 'FontSize', 14);
xlim([-60 10]);
ylim([-30 30]); caxis([-1e-7,1e-7]);
xticks([-60 -30 0]);
xticklabels({'60°W','30°W','0°'});
yticks([flip(-1*[0 10 20 30]) [10 20 30]]);
yticklabels({'30°S','20°S','10°S','0°','10°N','20°N','30°N'});
ax = gca; set(ax, 'FontSize', 15, 'Color', 0.6*[1 1 1], 'TickDir', 'out', 'DataAspectRatioMode', 'auto', 'PlotBoxAspectRatioMode', 'auto');
title('d)', 'Units', 'normalized', 'Position', [0, 1.02], 'HorizontalAlignment', 'left', 'FontSize', 15);

% -------------------------------------------------------------------------
% Subplot 5 (Row 2, Col 2 - Wide): Atlantic Heatfluxes Binned Plot
% -------------------------------------------------------------------------
ss(5) = subplot('Position', [px(2) py(2) axw2 axh]);

contourf(1:1461,latitude_a,adv_regiona,-1e-7:1e-9:1e-7,'edgecolor','none');
hold on;
contour(1:1461,latitude_a,sst_regiona, [-1e-2, 0], 'LineColor', 'black', 'LineStyle', '--');
contour(1:1461,latitude_a,sst_regiona, 0:100:100, 'LineColor', 'black', 'LineWidth', 1.2);
contour(1:1461,latitude_a,sst_regiona, [0 1e-2], 'LineColor', 'black');
cm = redblue(101);
colormap(cm);
cb = colorbar();
cb.Ticks = -4e-8:2e-8:4e-8; set(cb, 'AxisLocation', 'out', 'FontSize', 14);
ylim([-30 30]);
caxis([-5e-8,5e-8])
yticks([flip(-1*[0 10 20 30]) [10 20 30]]);
yticklabels({'30°S','20°S','10°S','0°','10°N','20°N','30°N'});
xticks([44,44+365,44+365*2,45+365*3]);
xticklabels({'2006','2007','2008','2009'});
ax = gca; set(ax, 'FontSize', 15);

title('e)', 'Units', 'normalized', 'Position', [0, 1.02], 'HorizontalAlignment', 'left', 'FontSize', 15);

% -------------------------------------------------------------------------
% Subplot 6 (Row 2, Col 3): Pacific Virtual Pot Temp Grad Binned Plot
% -------------------------------------------------------------------------
m = size(sst_regiona, 1);
n = size(sst_regiona,2);
r = zeros(m, 1);
p = zeros(m, 1);

for i = 1:m
    [R, P] = corrcoef(sst_regiona(i,:), adv_regiona(i,:));
    r(i) = R(1, 2);
    p(i) = P(1, 2);
end

% --- 底部 x 轴（相关系数）---
ax1 = axes('Position', [px(3) py(2) axw3 axh], 'XAxisLocation', 'bottom', 'YAxisLocation', 'left');
plot(ax1, movmean(r,12), latitude_a, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 4);
xlabel(ax1, 'Correlation coefficient', 'FontSize', 15);
ylim([-30 30]);
yticks([flip(-1*[0 10 20 30]) [10 20 30]]);
yticklabels({'30°S','20°S','10°S','0°','10°N','20°N','30°N'});
ax1.XColor = 'b';
ax1.YColor = 'k';
grid(ax1, 'on');
xlim(ax1, [-1.1, 1.1]);
set(ax1, 'FontSize', 15);

% --- 顶部 x 轴（p 值）---
ax2 = axes('Position', [px(3) py(2) axw3 axh], 'XAxisLocation', 'top', 'YAxisLocation', 'left', ...
           'Color', 'none', 'Box', 'off');
ax2.YAxis.Visible = 'off';
ax2.YLim = ax1.YLim;
hold(ax2, 'on');
plot(ax2, p, latitude_a, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 4);
%xlabel(ax2, 'p-value', 'FontSize', 15);
ax2.XColor = 'r';
xlim(ax2, [-0.15, 2.05]);
set(ax2, 'FontSize', 15);

xline(ax2, 0.05, 'k--', 'LineWidth', 1.5);
text(ax2, 0.05, double(max(latitude_a)*0.9), ' p = 0.05', 'Color', 'k', ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 15);
title('f)', 'Units', 'normalized', 'Position', [0, 1.02], 'HorizontalAlignment', 'left', 'FontSize', 15);

% Save figure
saveas(gcf, 'sst_total_heat_flux_2x3.fig');
close;