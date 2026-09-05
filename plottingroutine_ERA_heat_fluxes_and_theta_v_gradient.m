load era5_theta_v_height_coords.mat
load era5_theta_v_surface2m.mat

lon = ERA5_ThetaV.longitude;
lonp = lon+360;
latp = ERA5_ThetaV.latitude;

z_target = [2, 100:100:3000];
[nlon, nlat, nlevel, ntime] = size( ERA5_ThetaV.theta_v);
theta_v_2m = reshape(theta_v_sfc, nlon, nlat, 1, ntime);
clear theta_v_sfc;
theta_v = cat(3, theta_v_2m, ERA5_ThetaV.theta_v(:,:,2:end,:));
clear theta_v_2m ERA5_ThetaV
% Define height levels for gradient calculation (e.g., 500m to 1500m)
z_lower = 2;  % [m]
z_upper = 2000; % [m]

% Calculate vertical gradient between 500m and 1500m
[dth_dz, z_mid, stability] = calc_theta_v_gradient(theta_v, z_target, z_lower, z_upper, 3);
clear theta_v z_mid stability
dth_dz_mean_p = squeeze(nanmean(dth_dz(:,:,:,11:59),4));
clear dth_dz;
[I,J]=size(dth_dz_mean_p);
load etopo_globe.mat
landmaskp =size(dth_dz_mean_p);
for i=1:I
    for j=1:J
        lonind = find(lon>=lonp(i));  
        latind = find(lat<=latp(j));
        if topo(latind(1),lonind(1))<=0
            landmaskp(i,j)=topo(latind(1),lonind(1));
        else
            landmaskp(i,j)=NaN;
        end
    end
end    
clear topo

for i =1:I
    nanlocy =  find(landmaskp(i,:)~=landmaskp(i,:));
    dth_dz_mean_p(i,nanlocy,:)=NaN;
end

[I,J]=size(dth_dz_mean_p);
A = find(lonp>260);
B = find(latp>0);
for j=1:length(B)
    nanloc = find(dth_dz_mean_p(A,B(j))~=dth_dz_mean_p(A,B(j)));
    dth_dz_mean_p(A(nanloc(1)):I,B(j))=NaN;
end

dth_dz_mean_p = permute(dth_dz_mean_p,[2,1]);

dth_dz_mean_p1=movmean(dth_dz_mean_p,4,1,'omitnan');
dth_dz_mean_p=movmean(dth_dz_mean_p1,4,2,'omitnan');

dth_dz_mean_lowpassed1=movmean(dth_dz_mean_p,16,1,'omitnan');
dth_dz_mean_lowpassed2=movmean(dth_dz_mean_lowpassed1,16,2,'omitnan');clear dth_dz_mean_lowpassed1
dth_dz_mean_lowpassed3=movmean(dth_dz_mean_p-dth_dz_mean_lowpassed2,16,1,'omitnan');
dth_dz_mean_lowpassed4=movmean(dth_dz_mean_lowpassed3,16,2,'omitnan');clear dth_dz_mean_lowpassed3
dth_dz_meanp_highpassed=dth_dz_mean_p-(dth_dz_mean_lowpassed4+dth_dz_mean_lowpassed2);clear dth_dz_mean_lowpassed2 dth_dz_mean_lowpassed4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = 'era5_monthly_heatfluxes_2005-2009.nc';

slhf = ncread(filename, 'slhf')./86400;
sshf = ncread(filename, 'sshf')./86400;
ssr = ncread(filename, 'ssr')./86400;
str = ncread(filename, 'str')./86400;
slhf_mean = squeeze(nanmean(slhf(:,:,11:59), 3));clear slhf;
sshf_mean = squeeze(nanmean(sshf(:,:,11:59), 3));clear sshf;
ssr_mean = squeeze(nanmean(ssr(:,:,11:59), 3));clear ssr;
str_mean = squeeze(nanmean(str(:,:,11:59), 3));clear str;
totalhf_mean_p = slhf_mean+sshf_mean+ssr_mean+str_mean;clear slhf_mean sshf_mean ssr_mean str_mean

[I,J]=size(totalhf_mean_p);

for i =1:I
    nanlocy =  find(landmaskp(i,:)~=landmaskp(i,:));
    totalhf_mean_p(i,nanlocy,:)=NaN;
end

A = find(lonp>260);
B = find(latp>0);
for j=1:length(B)
    nanloc = find(totalhf_mean_p(A,B(j))~=totalhf_mean_p(A,B(j)));
    totalhf_mean_p(A(nanloc(1)):I,B(j))=NaN;
end

totalhf_mean_p = permute(totalhf_mean_p,[2,1]);

totalhf_mean_p1=movmean(totalhf_mean_p,4,1,'omitnan');
totalhf_mean_p=movmean(totalhf_mean_p1,4,2,'omitnan');clear totalhf_mean_p1

totalhf_mean_lowpassed1=movmean(totalhf_mean_p,16,1,'omitnan');
totalhf_mean_lowpassed2=movmean(totalhf_mean_lowpassed1,16,2,'omitnan');clear totalhf_mean_lowpassed1
totalhf_mean_lowpassed3=movmean(totalhf_mean_p-totalhf_mean_lowpassed2,16,1,'omitnan');
totalhf_mean_lowpassed4=movmean(totalhf_mean_lowpassed3,16,2,'omitnan');clear totalhf_mean_lowpassed3
totalhf_meanp_highpassed=totalhf_mean_p-(totalhf_mean_lowpassed4+totalhf_mean_lowpassed2);clear totalhf_mean_lowpassed2 totalhf_mean_lowpassed4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('D:\DATA\Ifremer\IfremerdailyTAUandSSTAllPacific.mat');
clear taux tauy;
sst=sst-273.15;

sst1=movmean(sst,4,1,'omitnan');
sst=movmean(sst1,4,2,'omitnan');
clear sst1

[I,J,K]=size(sst);
A = find(lons>260);
B = find(lats>0);
for j=1:length(B)
    nanloc = find(sst(A,B(j),1)~=sst(A,B(j),1));
    sst(A(nanloc(1)):I,B(j),:)=NaN;
end

load etopo_globe.mat
landmaskp =size(sst(:,:,1));
for i=1:I
    for j=1:J
        lonind = find(lon>=lons(i));  
        latind = find(lat<=lats(j));
        if topo(latind(1),lonind(1))<=0
            landmaskp(i,j)=topo(latind(1),lonind(1));
        else
            landmaskp(i,j)=NaN;
        end
    end
end    
clear topo

for i =1:I
    nanlocy =  find(landmaskp(i,:)~=landmaskp(i,:));
    sst(i,nanlocy,:)=NaN;
end

sst = permute(sst,[2,1,3]);

sstmean=nanmean(sst,3);clear sst
sstmeans_lowpassed1=movmean(sstmean,16,1,'omitnan');
sstmeans_lowpassed2=movmean(sstmeans_lowpassed1,16,2,'omitnan');clear sstmeans_lowpassed1
sstmeans_lowpassed3=movmean(sstmean-sstmeans_lowpassed2,16,1,'omitnan');
sstmeans_lowpassed4=movmean(sstmeans_lowpassed3,16,2,'omitnan');clear sstmeans_lowpassed3
sstmeans_highpassed=sstmean-(sstmeans_lowpassed4+sstmeans_lowpassed2);

%% ATLANTIC
load ATLera5_theta_v_height_coords.mat
load ATLera5_theta_v_surface2m.mat

lon = ERA5_ThetaV.longitude;
lona = lon+360;
lata = ERA5_ThetaV.latitude;

z_target = [2, 100:100:3000];
[nlon, nlat, nlevel, ntime] = size( ERA5_ThetaV.theta_v);
theta_v_2m = reshape(theta_v_sfc, nlon, nlat, 1, ntime);
clear theta_v_sfc;
theta_v = cat(3, theta_v_2m, ERA5_ThetaV.theta_v(:,:,2:end,:));
clear theta_v_2m ERA5_ThetaV
% Define height levels for gradient calculation (e.g., 500m to 1500m)
z_lower = 2;  % [m]
z_upper = 2000; % [m]

% Calculate vertical gradient between 500m and 1500m
[dth_dz, z_mid, stability] = calc_theta_v_gradient(theta_v, z_target, z_lower, z_upper, 3);
clear theta_v;
dth_dz_mean_a = squeeze(nanmean(dth_dz(:,:,:,11:59),4));
clear dth_dz;

dth_dz_mean_a = permute(dth_dz_mean_a(81:end,:,:),[2,1]);

dth_dz_mean_a1=movmean(dth_dz_mean_a,4,1,'omitnan');
dth_dz_mean_a=movmean(dth_dz_mean_a1,4,2,'omitnan');

dth_dz_mean_lowpassed1=movmean(dth_dz_mean_a,16,1,'omitnan');
dth_dz_mean_lowpassed2=movmean(dth_dz_mean_lowpassed1,16,2,'omitnan');clear dth_dz_mean_lowpassed1
dth_dz_mean_lowpassed3=movmean(dth_dz_mean_a-dth_dz_mean_lowpassed2,16,1,'omitnan');
dth_dz_mean_lowpassed4=movmean(dth_dz_mean_lowpassed3,16,2,'omitnan');clear dth_dz_mean_lowpassed3
dth_dz_meana_highpassed=dth_dz_mean_a-(dth_dz_mean_lowpassed4+dth_dz_mean_lowpassed2);clear dth_dz_mean_lowpassed2 dth_dz_mean_lowpassed4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = 'ATLera5_monthly_heatfluxes_2005-2009.nc';

slhf = ncread(filename, 'slhf')./86400;
sshf = ncread(filename, 'sshf')./86400;
ssr = ncread(filename, 'ssr')./86400;
str = ncread(filename, 'str')./86400;
slhf_mean = squeeze(nanmean(slhf(:,:,11:59), 3));clear slhf;
sshf_mean = squeeze(nanmean(sshf(:,:,11:59), 3));clear sshf;
ssr_mean = squeeze(nanmean(ssr(:,:,11:59), 3));clear ssr;
str_mean = squeeze(nanmean(str(:,:,11:59), 3));clear str;
totalhf_mean_a = slhf_mean+sshf_mean+ssr_mean+str_mean;clear slhf_mean sshf_mean ssr_mean str_mean

totalhf_mean_a = permute(totalhf_mean_a(81:end,:,:),[2,1]);

totalhf_mean_a1=movmean(totalhf_mean_a,4,1,'omitnan');
totalhf_mean_a=movmean(totalhf_mean_a1,4,2,'omitnan');clear totalhf_mean_a1

totalhf_mean_lowpassed1=movmean(totalhf_mean_a,16,1,'omitnan');
totalhf_mean_lowpassed2=movmean(totalhf_mean_lowpassed1,16,2,'omitnan');clear totalhf_mean_lowpassed1
totalhf_mean_lowpassed3=movmean(totalhf_mean_a-totalhf_mean_lowpassed2,16,1,'omitnan');
totalhf_mean_lowpassed4=movmean(totalhf_mean_lowpassed3,16,2,'omitnan');clear totalhf_mean_lowpassed3
totalhf_meana_highpassed=totalhf_mean_a-(totalhf_mean_lowpassed4+totalhf_mean_lowpassed2);clear totalhf_mean_lowpassed2 totalhf_mean_lowpassed4

lona = lona(81:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('D:\DATA\Ifremer\IfremerdailyWindstressAndSSTAtlantic.mat');
clear taux tauy;
sst=sst-273.15;

sst1=movmean(sst,4,1,'omitnan');
sst=movmean(sst1,4,2,'omitnan');
clear sst1

sst = permute(sst,[2,1,3]);

sstmeana=nanmean(sst,3);
clear sst;
sstmeans_lowpassed1=movmean(sstmeana,16,1,'omitnan');
sstmeans_lowpassed2=movmean(sstmeans_lowpassed1,16,2,'omitnan');clear sstmeans_lowpassed1
sstmeans_lowpassed3=movmean(sstmeana-sstmeans_lowpassed2,16,1,'omitnan');
sstmeans_lowpassed4=movmean(sstmeans_lowpassed3,16,2,'omitnan');clear sstmeans_lowpassed3
sstmeana_highpassed=sstmeana-(sstmeans_lowpassed4+sstmeans_lowpassed2);

dth_dz_meana_highpassed(isnan(flip(sstmeana,1)))=NaN;
totalhf_meana_highpassed(isnan(flip(sstmeana,1)))=NaN;
lona = lona-360;
%%
% --- MAIN PLOTTING SCRIPT ---
% =========================================================================
% Layout Parameters Setup
% =========================================================================
marg_h = [0.07, 0.04]; % Reduced bottom margin [bottom, top]
marg_w = [0.04, 0.04]; % Reduced left margin [left, right]
Nh = 2; 
Nw = 4;
gap = [0.07, 0.08, 0.07]; % [Reduced vert gap, Increased horiz gap 1, Increased horiz gap 2]

% Calculate heights and widths with safe proportions
axh  = (1 - sum(marg_h) - gap(1)) * 0.66; % Top panel height
axh1 = (1 - sum(marg_h) - gap(1)) * 0.33; % Bottom panel height
axw  = (1 - sum(marg_w) - 2*gap(2) - gap(3)) / Nw;

% Exact X and Y anchor positions
px = [marg_w(1), ...
      marg_w(1) + axw + gap(2), ...
      marg_w(1) + 2*axw + 2*gap(2), ...
      marg_w(1) + 3*axw + 2*gap(2) + gap(3)];

py = [1 - marg_h(2) - axh, ...
      marg_h(1)]; 

% -------------------------------------------------------------------------
% Subplot 1: Pacific heatfluxes & SST Contour (Panel a)
% -------------------------------------------------------------------------
ss(1) = subplot('Position', [px(1) py(1) axw axh]);
contourf(lonp, latp, flip(sstmeans_highpassed,1), -5e-1:5e-3:5e-1, 'LineStyle', 'None');
hold on;
contour(lonp, latp,totalhf_meanp_highpassed,0:100:100,'LineColor','black','LineWidth',1.2);
contour(lonp, latp,totalhf_meanp_highpassed,[-1e-0,0] ,'LineColor','black','LineStyle','--');
contour(lonp, latp,totalhf_meanp_highpassed,0:1e-0:1e-0,'LineColor','black');
cm = redblue(101);
colormap(gca, cm);
worldmap3(2);
cb = colorbar();
cb.Ticks = -0.08:0.04:0.08;
set(cb, 'AxisLocation', 'out', 'FontSize', 14); % Increased colorbar font size
xlim([210 290]);
ylim([-30 30]);
caxis([-10e-2, 10e-2]);
xticks([210 240 270]);
xticklabels({'150¡ãW','120¡ãW','90¡ãW'});
yticks([flip(-1*[0 10 20 30]) [10 20 30]]);
yticklabels({'30¡ãS','20¡ãS','10¡ãS','0¡ã','10¡ãN','20¡ãN','30¡ãN'});
ax = gca;
set(ax, 'FontSize', 14, 'Color', 0.6*[1 1 1], 'TickDir', 'out', ...
        'DataAspectRatioMode', 'auto', 'PlotBoxAspectRatioMode', 'auto', ...
        'Position', [px(1) py(1) axw axh]);
ax.YAxis.TickLabelGapOffset = 0;
ax.XAxis.TickLabelGapOffset = 0;
title('a)', 'Units', 'normalized', 'Position', [0, 1.02], 'HorizontalAlignment', 'left', 'FontSize', 15);

% -------------------------------------------------------------------------
% Subplot 2: Pacific Sheatfluxes & SST Binned Plot (Panel e)
% -------------------------------------------------------------------------
ss(2) = subplot('Position', [px(1) py(2) axw axh1]);
x_min = -10; x_max = 10;
[bx, by, bstd, slope, x_fit, y_fit] = compute_binned_stats( ...
    totalhf_meanp_highpassed .* 1e0, flip(sstmeans_highpassed,1) .* 1e1,  x_min, x_max, 20, 4);

line([x_min*2, x_max*2], [0, 0], 'Color', 'k', 'LineWidth', 0.8); hold on;
plot(x_fit, y_fit, 'k-', 'LineWidth', 1.2);
errorbar(bx, by, bstd, 'o', 'Color', 'k', 'MarkerSize', 4, ...
    'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'LineWidth', 1.0, 'CapSize', 3);

xlim([x_min*1.5, x_max*1.5]);
ylim([-3,3]);
yticks(-2:1:2);
xlabel('total heat flux (W m$^{-2}$)', 'Interpreter', 'latex');
yl = ylabel('SST ($10^{-1}$ $^\circ\mathrm{C}$ )', 'Interpreter', 'latex');
yl.Units = 'normalized';
yl.Position(1) = -0.08;
ax = gca;
set(ax, 'FontSize', 14, 'Box', 'on', 'TickDir', 'in', 'XMinorTick', 'on', 'YMinorTick', 'on');
ax.YAxis.TickLabelGapOffset = 0; 
text(0.65, 0.88, sprintf('s = %.2f', slope), 'Units', 'normalized', 'FontSize', 14, 'FontName', 'Helvetica');
title('e)', 'Units', 'normalized', 'Position', [0, 1.02], 'HorizontalAlignment', 'left', 'FontSize', 15);

% -------------------------------------------------------------------------
% Subplot 3: Atlantic heatfluxes & SST Contour (Panel b)
% -------------------------------------------------------------------------
ss(5) = subplot('Position', [px(2) py(1) axw axh]);
contourf(lona, lata, flip(sstmeana_highpassed,1), -5e-1:5e-3:5e-1, 'LineStyle', 'None');
hold on;
contour(lona, lata,totalhf_meana_highpassed,0:100:100,'LineColor','black','LineWidth',1.2);
contour(lona, lata,totalhf_meana_highpassed,[-1e-0,0] ,'LineColor','black','LineStyle','--');
contour(lona, lata,totalhf_meana_highpassed,0:1e-0:1e-0,'LineColor','black');
cm = redblue(101);
colormap(gca, cm);
worldmap3Atl(2);
cb = colorbar();
cb.Ticks = -0.08:0.04:0.08;
set(cb, 'AxisLocation', 'out', 'FontSize', 14); % Increased colorbar font size
xlim([-60 10]);
ylim([-30 30]);
caxis([-10e-2, 10e-2]);
xticks([-60 -30 0]);
xticklabels({'60¡ãW','30¡ãW','0¡ã'});
yticks([flip(-1*[0 10 20 30]) [10 20 30]]);
yticklabels({'30¡ãS','20¡ãS','10¡ãS','0¡ã','10¡ãN','20¡ãN','30¡ãN'});
title('b)', 'Units', 'normalized', 'Position', [0, 1.02], 'HorizontalAlignment', 'left', 'FontSize', 15);
ax = gca;
set(ax, 'FontSize', 14, 'Color', 0.6*[1 1 1], 'TickDir', 'out', ...
        'DataAspectRatioMode', 'auto', 'PlotBoxAspectRatioMode', 'auto', ...
        'Position', [px(2) py(1) axw axh]);
ax.YAxis.TickLabelGapOffset = 0;
ax.XAxis.TickLabelGapOffset = 0;

% -------------------------------------------------------------------------
% Subplot 4: Atlantic heatfluxes & SST Binned Plot (Panel f)
% -------------------------------------------------------------------------
ss(6) = subplot('Position', [px(2) py(2) axw axh1]);
x_min = -9; x_max = 9;
[bx, by, bstd, slope, x_fit, y_fit] = compute_binned_stats( ...
    totalhf_meana_highpassed .* 1e0, flip(sstmeana_highpassed,1) .* 1e1, x_min, x_max, 20,5);

line([x_min*2, x_max*2], [0, 0], 'Color', 'k', 'LineWidth', 0.8); hold on;
plot(x_fit, y_fit, 'k-', 'LineWidth', 1.2);
errorbar(bx, by, bstd, 'o', 'Color', 'k', 'MarkerSize', 4, ...
    'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'LineWidth', 1.0, 'CapSize', 3);

xlim([x_min*15/9, x_max*15/9]);
ylim([-3,3]);
yticks(-2:1:2);
xlabel('total heat flux (W m$^{-2}$)', 'Interpreter', 'latex');
yl = ylabel('SST ($10^{-1}$ $^\circ\mathrm{C}$ )', 'Interpreter', 'latex');
ax = gca;
set(ax, 'FontSize', 14, 'Box', 'on', 'TickDir', 'in', 'XMinorTick', 'on', 'YMinorTick', 'on');
ax.YAxis.TickLabelGapOffset = 0;
text(0.65, 0.88, sprintf('s = %.2f', slope), 'Units', 'normalized', 'FontSize', 14, 'FontName', 'Helvetica');
title('f)', 'Units', 'normalized', 'Position', [0, 1.02], 'HorizontalAlignment', 'left', 'FontSize', 15);

% -------------------------------------------------------------------------
% Subplot 5: Pacific virtual pot temp grad & SST Contour (Panel c)
% -------------------------------------------------------------------------
ss(3) = subplot('Position', [px(3) py(1) axw axh]);
contourf(lonp, latp, totalhf_meanp_highpassed, -25e-0:5e-1:25e-0, 'LineStyle', 'None');
hold on;
contour(lonp, latp,dth_dz_meanp_highpassed,0:100:100,'LineColor','black','LineWidth',1.2);
contour(lonp, latp,dth_dz_meanp_highpassed,[-5e-6,0] ,'LineColor','black','LineStyle','--');
contour(lonp, latp,dth_dz_meanp_highpassed,0:5e-6:5e-6,'LineColor','black');
cm = redblue(101);
colormap(gca, cm);
worldmap3(2);
cb = colorbar();
cb.Ticks = -4:2:4;
set(cb, 'AxisLocation', 'out', 'FontSize', 14); % Increased colorbar font size
xlim([210 290]);
ylim([-30 30]);
caxis([-5e-0, 5e-0]);
xticks([210 240 270]);
xticklabels({'150¡ãW','120¡ãW','90¡ãW'});
yticks([flip(-1*[0 10 20 30]) [10 20 30]]);
yticklabels({'30¡ãS','20¡ãS','10¡ãS','0¡ã','10¡ãN','20¡ãN','30¡ãN'});
title('c)', 'Units', 'normalized', 'Position', [0, 1.02], 'HorizontalAlignment', 'left', 'FontSize', 15);
ax = gca;
set(ax, 'FontSize', 14, 'Color', 0.6*[1 1 1], 'TickDir', 'out', ...
        'DataAspectRatioMode', 'auto', 'PlotBoxAspectRatioMode', 'auto', ...
        'Position', [px(3) py(1) axw axh]);
ax.YAxis.TickLabelGapOffset = 0;
ax.XAxis.TickLabelGapOffset = 0;

% -------------------------------------------------------------------------
% Subplot 6: Pacific virtual pot temp grad & SST Binned Plot (Panel g)
% -------------------------------------------------------------------------
ss(4) = subplot('Position', [px(3) py(2) axw axh1]);

x_min = -30; x_max = 30;
[bx, by, bstd, slope, x_fit, y_fit] = compute_binned_stats( ...
    dth_dz_meanp_highpassed .* 1e6,totalhf_meanp_highpassed .* 1e0,  x_min, x_max, 20, 4);

line([x_min*2, x_max*2], [0, 0], 'Color', 'k', 'LineWidth', 0.8); hold on;
plot(x_fit, y_fit, 'k-', 'LineWidth', 1.2);
errorbar(bx, by, bstd, 'o', 'Color', 'k', 'MarkerSize', 4, ...
    'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'LineWidth', 1.0, 'CapSize', 3);

xlim([x_min*1.5, x_max*1.5]);
ylim([-10,10]);
yticks(-10:5:10);
xlabel('${\partial \theta_v} / {\partial z}$ ($10^{-6}$ $^\circ\mathrm{C}$ m$^{-1}$)', 'Interpreter', 'latex');
yl = ylabel('total heat flux (W m$^{-2}$)', 'Interpreter', 'latex');
yl.Units = 'normalized';
ax = gca;
set(ax, 'FontSize', 14, 'Box', 'on', 'TickDir', 'in', 'XMinorTick', 'on', 'YMinorTick', 'on');
ax.YAxis.TickLabelGapOffset = 0;
text(0.05, 0.88, sprintf('s = %.2f', slope), 'Units', 'normalized', 'FontSize', 14, 'FontName', 'Helvetica');
title('g)', 'Units', 'normalized', 'Position', [0, 1.02], 'HorizontalAlignment', 'left', 'FontSize', 15);

% -------------------------------------------------------------------------
% Subplot 7: Atlantic virtual pot temp grad & SST Contour (Panel d)
% -------------------------------------------------------------------------
ss(7) = subplot('Position', [px(4) py(1) axw axh]);
contourf(lona, lata,totalhf_meana_highpassed, -25e-0:5e-1:25e-0, 'LineStyle', 'None');
hold on;
contour(lona, lata, dth_dz_meana_highpassed, [-5e-6,0], 'LineColor', 'black', 'LineStyle', '--');
contour(lona, lata, dth_dz_meana_highpassed, 0:100:100, 'LineColor', 'black', 'LineWidth', 1.2);
contour(lona, lata, dth_dz_meana_highpassed, 0:5e-6:5e-6, 'LineColor', 'black');
cm = redblue(101);
colormap(gca, cm);
worldmap3Atl(2);
cb = colorbar();
cb.Ticks = -4:2:4;
set(cb, 'AxisLocation', 'out', 'FontSize', 14); % Increased colorbar font size
xlim([-60 10]);
ylim([-30 30]);
caxis([-5e-0, 5e-0]);
xticks([-60 -30 0]);
xticklabels({'60¡ãW','30¡ãW','0¡ã'});
yticks([flip(-1*[0 10 20 30]) [10 20 30]]);
yticklabels({'30¡ãS','20¡ãS','10¡ãS','0¡ã','10¡ãN','20¡ãN','30¡ãN'});
title('d)', 'Units', 'normalized', 'Position', [0, 1.02], 'HorizontalAlignment', 'left', 'FontSize', 15);
ax = gca;
set(ax, 'FontSize', 14, 'Color', 0.6*[1 1 1], 'TickDir', 'out', ...
        'DataAspectRatioMode', 'auto', 'PlotBoxAspectRatioMode', 'auto', ...
        'Position', [px(4) py(1) axw axh]);
ax.YAxis.TickLabelGapOffset = 0;
ax.XAxis.TickLabelGapOffset = 0;

% -------------------------------------------------------------------------
% Subplot 8: Atlantic virtual pot temp grad & SST Binned Plot (Panel h)
% -------------------------------------------------------------------------
ss(8) = subplot('Position', [px(4) py(2) axw axh1]);
x_min = -30; x_max = 30;
[bx, by, bstd, slope, x_fit, y_fit] = compute_binned_stats( ...
     dth_dz_meana_highpassed .* 1e6, totalhf_meana_highpassed,x_min, x_max, 20, 4);

line([x_min*2, x_max*2], [0, 0], 'Color', 'k', 'LineWidth', 0.8); hold on;
plot(x_fit, y_fit, 'k-', 'LineWidth', 1.2);
errorbar(bx, by, bstd, 'o', 'Color', 'k', 'MarkerSize', 4, ...
    'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'LineWidth', 1.0, 'CapSize', 3);

xlim([x_min*1.5, x_max*1.5]);
ylim([-10,10]);
yticks(-10:5:10);
xlabel('${\partial \theta_v} / {\partial z}$ ($10^{-6}$ $^\circ\mathrm{C}$ m$^{-1}$)', 'Interpreter', 'latex');
yl = ylabel('total heat flux (W m$^{-2}$)', 'Interpreter', 'latex');
ax = gca;
set(ax, 'FontSize', 14, 'Box', 'on', 'TickDir', 'in', 'XMinorTick', 'on', 'YMinorTick', 'on');
ax.YAxis.TickLabelGapOffset = 0;
text(0.05, 0.88, sprintf('s = %.2f', slope), 'Units', 'normalized', 'FontSize', 14, 'FontName', 'Helvetica');
title('h)', 'Units', 'normalized', 'Position', [0, 1.02], 'HorizontalAlignment', 'left', 'FontSize', 15);
saveas(gcf,'sst total_heat_flux and thetav grad.png');
close;
