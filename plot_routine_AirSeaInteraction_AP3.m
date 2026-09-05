load('IfremerdailyTAUandSSTAllPacific.mat');

taux(isnan(sst))=NaN;
tauy(isnan(sst))=NaN;
sst=sst-273.15;

taux1=movmean(taux,4,1,'omitnan');
taux=movmean(taux1,4,2,'omitnan');
tauy1=movmean(tauy,4,1,'omitnan');
tauy=movmean(tauy1,4,2,'omitnan');
sst1=movmean(sst,4,1,'omitnan');
sst=movmean(sst1,4,2,'omitnan');
clear taux1 tauy1 sst1

[I,J,K]=size(taux);
A = find(lons>260);
B = find(lats>0);
for j=1:length(B)
    nanloc = find(taux(A,B(j),1)~=taux(A,B(j),1));
    taux(A(nanloc(1)):I,B(j),:)=NaN;
    tauy(A(nanloc(1)):I,B(j),:)=NaN;
    sst(A(nanloc(1)):I,B(j),:)=NaN;
end

load etopo_globe.mat
landmaskp =size(taux(:,:,1));
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
    taux(i,nanlocy,:)=NaN;
    tauy(i,nanlocy,:)=NaN;
    sst(i,nanlocy,:)=NaN;
end


taux = permute(taux,[2,1,3]);
tauy = permute(tauy,[2,1,3]);
sst = permute(sst,[2,1,3]);
tauxmean=nanmean(taux,3);
tauymean=nanmean(tauy,3);

% partial derivatives of taux and tauy
[txx,txy]=zh_grad2(tauxmean,lons,lats);
clear txx 
[tyx,tyy]=zh_grad2(tauymean,lons,lats);
clear tyy 
% curl of tau
wsc=tyx-txy;
clear txy tyx

sstmean=nanmean(sst,3);
wscmeans=nanmean(wsc,3);
% gradient of sst
[sstx,ssty]=zh_grad2(sst,lons,lats);

[sstxx,sstxy]=zh_grad2(sstx,lons,lats);
clear sstxy;
[sstyx,sstyy]=zh_grad2(ssty,lons,lats);
clear sstyx;
% SST laplacian
sstlap=sstxx+sstyy;clear sstxx sstyy
sstlapmean_p=nanmean(sstlap,3);clear sstlap

sstlapmean_lowpassed1=movmean(sstlapmean_p,16,1,'omitnan');
sstlapmean_lowpassed2=movmean(sstlapmean_lowpassed1,16,2,'omitnan');clear sstlapmean_lowpassed1
sstlapmean_lowpassed3=movmean(sstlapmean_p-sstlapmean_lowpassed2,16,1,'omitnan');
sstlapmean_lowpassed4=movmean(sstlapmean_lowpassed3,16,2,'omitnan');clear sstlapmean_lowpassed3
sstlapmeanp_highpassed=sstlapmean_p-(sstlapmean_lowpassed4+sstlapmean_lowpassed2);clear sstlapmean_lowpassed2 sstlapmean_lowpassed4

taudir=atan2(tauy,taux);
taumag=sqrt(taux.^2+tauy.^2);clear taux tauy
% sst gradient direction and magnitude
sstgraddir=atan2(ssty,sstx);
sstgradmag=sqrt(sstx.^2+ssty.^2);clear sstx ssty
%crosswind sst gradient
crosswindsstgrad = sstgradmag.*sin(taudir-sstgraddir);
cwsgs_mean=nanmean(crosswindsstgrad,3);clear crosswindsstgrad

wscmeans_lowpassed1=movmean(wscmeans,16,1,'omitnan');
wscmeans_lowpassed2=movmean(wscmeans_lowpassed1,16,2,'omitnan');clear wscmeans_lowpassed1
wscmeans_lowpassed3=movmean(wscmeans-wscmeans_lowpassed2,16,1,'omitnan');
wscmeans_lowpassed4=movmean(wscmeans_lowpassed3,16,2,'omitnan');clear wscmeans_lowpassed3
wscmeans_highpassed=wscmeans-(wscmeans_lowpassed4+wscmeans_lowpassed2);


sstmeans_lowpassed1=movmean(sstmean,16,1,'omitnan');
sstmeans_lowpassed2=movmean(sstmeans_lowpassed1,16,2,'omitnan');clear sstmeans_lowpassed1
sstmeans_lowpassed3=movmean(sstmean-sstmeans_lowpassed2,16,1,'omitnan');
sstmeans_lowpassed4=movmean(sstmeans_lowpassed3,16,2,'omitnan');clear sstmeans_lowpassed3
sstmeans_highpassed=sstmean-(sstmeans_lowpassed4+sstmeans_lowpassed2);



cwsgsmean_lowpassed1=movmean(cwsgs_mean,16,1,'omitnan');
cwsgsmean_lowpassed2=movmean(cwsgsmean_lowpassed1,16,2,'omitnan');
cwsgsmean_lowpassed3=movmean(cwsgs_mean-cwsgsmean_lowpassed2,16,1,'omitnan');
cwsgsmean_lowpassed4=movmean(cwsgsmean_lowpassed3,16,2,'omitnan');
cwsgsmean_highpassed=cwsgs_mean-(cwsgsmean_lowpassed4+cwsgsmean_lowpassed2);

wscmeans_1d=reshape(wscmeans_highpassed,[],1);
cwsgsmean_1d=reshape(cwsgsmean_highpassed,[],1);
sstlapmeanp_1d=reshape(sstlapmeanp_highpassed,[],1);

[J,I,K]=size(cwsgsmean_highpassed);
A = find(lons>260);
B = find(lats>0);
for j=1:length(B)
    nanloc = find(cwsgsmean_highpassed(B(j),A)~=cwsgsmean_highpassed(B(j),A));
    cwsgsmean_highpassed(B(j),A(nanloc(1)):I)=0;
    wscmeans_highpassed(B(j),A(nanloc(1)):I)=0;
    sstlapmeanp_highpassed(B(j),A(nanloc(1)):I)=0;
end

for i =1:I
    nanlocy =  find(landmaskp(i,:)~=landmaskp(i,:));
    cwsgsmean_highpassed(nanlocy,i)=NaN;
    wscmeans_highpassed(nanlocy,i)=NaN;
    sstlapmeanp_highpassed(nanlocy,i)=NaN;
end
lonp = lons;
latp = lats;



%%
load('IfremerdailyWindstressAndSSTAtlantic.mat');
[I,J,K]=size(taux);
taux(isnan(sst))=NaN;
tauy(isnan(sst))=NaN;

lona = lons;clear lons
lata = lats;clear lats

taux = permute(taux,[2,1,3]);
tauy = permute(tauy,[2,1,3]);
sst = permute(sst,[2,1,3]);

taux1=movmean(taux,4,1,'omitnan');
taux=movmean(taux1,4,2,'omitnan');
tauy1=movmean(tauy,4,1,'omitnan');
tauy=movmean(tauy1,4,2,'omitnan');
sst1=movmean(sst,4,1,'omitnan');
sst=movmean(sst1,4,2,'omitnan');
clear taux1 tauy1 sst1
tauxmean=nanmean(taux,3);
tauymean=nanmean(tauy,3);

% partial derivatives of taux and tauy
[txx,txy]=zh_grad2(tauxmean,lona,lata);
clear txx 
[tyx,tyy]=zh_grad2(tauymean,lona,lata);
clear tyy 
% curl of tau
wsc=tyx-txy;
clear txy tyx

sst=sst-273.15;
sstmean=nanmean(sst,3);
wscmeana=nanmean(wsc,3);
% gradient of sst
[sstx,ssty]=zh_grad2(sst,lona,lata);

[sstxx,sstxy]=zh_grad2(sstx,lona,lata);
clear sstxy;
[sstyx,sstyy]=zh_grad2(ssty,lona,lata);
clear sstyx;
% SST laplacian
sstlap=sstxx+sstyy;clear sstxx sstyy
sstlapmean_a=nanmean(sstlap,3);clear sstlap

sstlapmean_lowpassed1=movmean(sstlapmean_a,16,1,'omitnan');
sstlapmean_lowpassed2=movmean(sstlapmean_lowpassed1,16,2,'omitnan');clear sstlapmean_lowpassed1
sstlapmean_lowpassed3=movmean(sstlapmean_a-sstlapmean_lowpassed2,16,1,'omitnan');
sstlapmean_lowpassed4=movmean(sstlapmean_lowpassed3,16,2,'omitnan');clear sstlapmean_lowpassed3
sstlapmeana_highpassed=sstlapmean_a-(sstlapmean_lowpassed4+sstlapmean_lowpassed2);clear sstlapmean_lowpassed2 sstlapmean_lowpassed4

taudir=atan2(tauy,taux);
taumag=sqrt(taux.^2+tauy.^2);clear taux tauy
% sst gradient direction and magnitude
sstgraddir=atan2(ssty,sstx);
sstgradmag=sqrt(sstx.^2+ssty.^2);clear sstx ssty
%crosswind sst gradient
crosswindsstgrad = sstgradmag.*sin(taudir-sstgraddir);clear taudir sstgraddir sstgradmag
cwsga_mean=nanmean(crosswindsstgrad,3);clear crosswindsstgrad

wscmeana_lowpassed1=movmean(wscmeana,16,1,'omitnan');
wscmeana_lowpassed2=movmean(wscmeana_lowpassed1,16,2,'omitnan');clear wscmeann_lowpassed1
wscmeana_lowpassed3=movmean(wscmeana-wscmeana_lowpassed2,16,1,'omitnan');
wscmeana_lowpassed4=movmean(wscmeana_lowpassed3,16,2,'omitnan');clear wscmeann_lowpassed3
wscmeana_highpassed=wscmeana-(wscmeana_lowpassed4+wscmeana_lowpassed2);clear wscmeana_lowpassed2 wscmeana_lowpassed4

cwsgamean_lowpassed1=movmean(cwsga_mean,16,1,'omitnan');
cwsgamean_lowpassed2=movmean(cwsgamean_lowpassed1,16,2,'omitnan');clear cwsgamean_lowpassed1
cwsgamean_lowpassed3=movmean(cwsga_mean-cwsgamean_lowpassed2,16,1,'omitnan');
cwsgamean_lowpassed4=movmean(cwsgamean_lowpassed3,16,2,'omitnan');clear cwsgamean_lowpassed3
cwsgamean_highpassed=cwsga_mean-(cwsgamean_lowpassed4+cwsgamean_lowpassed2);clear cwsgamean_lowpassed2 cwsgamean_lowpassed4

wscmeana_1d=reshape(wscmeana_highpassed,[],1);
cwsgamean_1d=reshape(cwsgamean_highpassed,[],1);
sstlapamean_1d=reshape(sstlapmeana_highpassed,[],1);
%plot crosswind sst gradient
ss(1)=subplot(1,2,1);
contourf(lona,lata,wscmeana_highpassed,-1e-7:1e-9:1e-7,'edgecolor','None');
hold on;
xlabel('Longitude');
ylabel('Latiutude');
cm = redblue(101);
colormap(cm);
colorbar();
worldmap3Atl(2);
xlim([-60 10]);
ylim([-30 30]);
title('wscmeana highpassed');

ss(2)=subplot(1,2,2);
contourf(lona,lata,cwsgamean_highpassed,-1e-5:1e-7:1e-5,'LineStyle','None');
hold on;
ylabel('Latitude');
xlabel('Longitude');
cm = redblue(101);
colormap(cm);
colorbar();
worldmap3Atl(2);
xlim([-60 10]);
ylim([-30 30]);
title('Mean Crosswind SST gradient highpassed');

% --- MAIN PLOTTING SCRIPT ---

marg_h = [0.06 0.04]; 
marg_w = [0.04 0.02];
Nh = 3;
Nw = 4;
gap = [0.055 0.035 0.045];

axh  = (1 - sum(marg_h) - 1*gap(1)) / 1.5; 
axw  = (1 - sum(marg_w) - 2*gap(2) - gap(3)) / 4;
axh1 = (1 - sum(marg_h) - 1*gap(1)) / Nh; 
axw1 = (1 - sum(marg_w) - 2*gap(2) - gap(3)) / 4;

px = [marg_w(1), ...
      marg_w(1) + axw + gap(2), ...
      marg_w(1) + axw + gap(2) + gap(3) + axw1, ...
      marg_w(1) + 2*axw + 2*gap(2) + gap(3) + axw1];
py = [1 - marg_h(2) - axh, ...
      1 - marg_h(2) - axh - gap(1) - axh1]; 

set(gcf, 'Units', 'centimeters', 'Position', [50, -5, 38.5, 25.5]);

% -------------------------------------------------------------------------
% Subplot 1: Pacific WSC & SST Laplacian Contour (Panel a)
% -------------------------------------------------------------------------
ss(1) = subplot('Position', [px(1) py(1) axw axh]);
contourf(lonp, latp, wscmeans_highpassed, -5e-7:0.5e-9:5e-7, 'LineStyle', 'None');
hold on;
contour(lonp, latp, sstlapmeanp_highpassed, 0:100:100, 'LineColor', 'black', 'LineWidth', 1.2);
cm = redblue(101);
colorbar();
colormap(gca, cm);
worldmap3(2);
xlim([210 290]);
ylim([-30 30]);
caxis([-3e-8, 3e-8]);
xticks([210 240 270]);
xticklabels({'150°„W','120°„W','90°„W'});
yticks([flip(-1*[0 10 20 30]) [10 20 30]]);
yticklabels({'30°„S','20°„S','10°„S','0°„','10°„N','20°„N','30°„N'});
ax = gca;
set(ax, 'FontSize', 15, 'Color', 0.6*[1 1 1], 'TickDir', 'out'); 
ax.YAxis.TickLabelGapOffset = -2.5;
ax.XAxis.TickLabelGapOffset = -4;
title('a)', 'Units', 'normalized', 'Position', [0, 1.005], 'HorizontalAlignment', 'left', 'FontSize', 15);

% Get actual rendered position of top panel 'a'
drawnow;
pos1 = ss(1).Position;

% -------------------------------------------------------------------------
% Subplot 2: Pacific SST Laplacian vs Curl Binned Plot (Panel e)
% -------------------------------------------------------------------------
ss(2) = subplot('Position', [pos1(1) py(2) pos1(3) axh1]);
x_min = -.4; x_max = .4;
[bx, by, bstd, slope, x_fit, y_fit] = compute_binned_stats( ...
    sstlapmeanp_highpassed .* 1e10, wscmeans_highpassed .* 1e7, x_min, x_max, 20, 4);

line([x_min*2, x_max*2], [0, 0], 'Color', 'k', 'LineWidth', 0.8); hold on;
plot(x_fit, y_fit, 'k-', 'LineWidth', 1.2);
errorbar(bx, by, bstd, 'o', 'Color', 'k', 'MarkerSize', 4, ...
    'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'LineWidth', 1.0, 'CapSize', 3);

xlim([x_min*1.5, x_max*1.5]);
ylim([-.5,.5]);
yticks(-0.4:0.2:0.4);
xlabel('$\nabla^{2} T$ (10$^{-10}$ $^\circ\mathrm{C}$ m$^{-2}$)', 'Interpreter', 'latex');
yl = ylabel('WSC (N m$^{-3} \times 10^7$)', 'Interpreter', 'latex');
yl.Units = 'normalized';
yl.Position(1) = -0.14; % Shift x-position closer to y-axis
ax = gca;
set(ax, 'FontSize', 13, 'Box', 'on', 'TickDir', 'in', 'XMinorTick', 'on', 'YMinorTick', 'on');
text(0.05, 0.88, sprintf('s = %.2f', slope), 'Units', 'normalized', 'FontSize', 13, 'FontName', 'Helvetica');
title('e)', 'Units', 'normalized', 'Position', [0, 1.005], 'HorizontalAlignment', 'left', 'FontSize', 15);

% -------------------------------------------------------------------------
% Subplot 3: Pacific WSC & Crosswind SST Grad Contour (Panel c)
% -------------------------------------------------------------------------
ss(3) = subplot('Position', [px(3) py(1) axw axh]);
contourf(lonp, latp, wscmeans_highpassed, -5e-7:1e-9:5e-7, 'LineStyle', 'None');
hold on;
contour(lonp, latp, cwsgsmean_highpassed, [-2e-7, 0], 'LineColor', 'black', 'LineStyle', '--');
contour(lonp, latp, cwsgsmean_highpassed, 0:100:100, 'LineColor', 'black', 'LineWidth', 1.2);
contour(lonp, latp, cwsgsmean_highpassed, 0:2e-7:2e-7, 'LineColor', 'black');
cm = redblue(101);
colorbar();
colormap(gca, cm);
worldmap3(2);
xlim([210 290]);
ylim([-30 30]);
caxis([-3e-8, 3e-8]);
xticks([210 240 270]);
xticklabels({'150°„W','120°„W','90°„W'});
yticks([flip(-1*[0 10 20 30]) [10 20 30]]);
yticklabels({'30°„S','20°„S','10°„S','0°„','10°„N','20°„N','30°„N'});
title('c)', 'Units', 'normalized', 'Position', [0, 1.005], 'HorizontalAlignment', 'left', 'FontSize', 11);
ax = gca;
set(ax, 'FontSize', 15, 'Color', 0.6*[1 1 1], 'TickDir', 'out');
ax.YAxis.TickLabelGapOffset = -2.5;
ax.XAxis.TickLabelGapOffset = -4;

% Get actual rendered position of top panel 'c'
drawnow;
pos3 = ss(3).Position;

% -------------------------------------------------------------------------
% Subplot 4: Pacific Crosswind SST Grad vs Curl Binned Plot (Panel g)
% -------------------------------------------------------------------------
ss(4) = subplot('Position', [pos3(1) py(2) pos3(3) axh1]);
x_min = -0.1; x_max = 0.1;
[bx, by, bstd, slope, x_fit, y_fit] = compute_binned_stats( ...
    cwsgsmean_highpassed .* 1e5, wscmeans_highpassed .* 1e7, x_min, x_max, 20, 4);

line([x_min*2, x_max*2], [0, 0], 'Color', 'k', 'LineWidth', 0.8); hold on;
plot(x_fit, y_fit, 'k-', 'LineWidth', 1.2);
errorbar(bx, by, bstd, 'o', 'Color', 'k', 'MarkerSize', 4, ...
    'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'LineWidth', 1.0, 'CapSize', 3);

xlim([x_min*1.5, x_max*1.5]);
ylim([-.35,.35]);
yticks(-0.3:0.1:0.3);
xlabel('crosswind $\nabla T$ ($10^{-5}$ $^\circ\mathrm{C}$ m$^{-1}$)', 'Interpreter', 'latex');
ylabel('WSC (N m$^{-3} \times 10^7$)', 'Interpreter', 'latex');
ax = gca;
set(ax, 'FontSize', 13, 'Box', 'on', 'TickDir', 'in', 'XMinorTick', 'on', 'YMinorTick', 'on');
text(0.05, 0.88, sprintf('s = %.2f', slope), 'Units', 'normalized', 'FontSize', 13, 'FontName', 'Helvetica');
title('g)', 'Units', 'normalized', 'Position', [0, 1.005], 'HorizontalAlignment', 'left', 'FontSize', 15);

% -------------------------------------------------------------------------
% Subplot 5: Atlantic WSC & SST Laplacian Contour (Panel b)
% -------------------------------------------------------------------------
ss(5) = subplot('Position', [px(2) py(1) axw axh]);
contourf(lona, lata, wscmeana_highpassed, -5e-7:1e-9:5e-7, 'LineStyle', 'None');
hold on;
contour(lona, lata, sstlapmeana_highpassed, 0:100:100, 'LineColor', 'black', 'LineWidth', 1.2);
cm = redblue(101);
colorbar();
colormap(gca, cm);
worldmap3Atl(2);
xlim([-60 10]);
ylim([-30 30]);
caxis([-3e-8, 3e-8]);
xticks([-60 -30 0]);
xticklabels({'60°„W','30°„W','0°„'});
yticks([flip(-1*[0 10 20 30]) [10 20 30]]);
yticklabels({'30°„S','20°„S','10°„S','0°„','10°„N','20°„N','30°„N'});
title('b)', 'Units', 'normalized', 'Position', [0, 1], 'HorizontalAlignment', 'left');
ax = gca;
set(ax, 'FontSize', 15, 'Color', 0.6*[1 1 1], 'TickDir', 'out');
ax.YAxis.TickLabelGapOffset = -2.5;
ax.XAxis.TickLabelGapOffset = -4;

% Get actual rendered position of top panel 'b'
drawnow;
pos5 = ss(5).Position;

% -------------------------------------------------------------------------
% Subplot 6: Atlantic SST Laplacian vs Curl Binned Plot (Panel f)
% -------------------------------------------------------------------------
ss(6) = subplot('Position', [pos5(1) py(2) pos5(3) axh1]);
x_min = -0.4; x_max = 0.4;
[bx, by, bstd, slope, x_fit, y_fit] = compute_binned_stats( ...
    sstlapmeana_highpassed .* 1e10, wscmeana_highpassed .* 1e7, x_min, x_max, 20, 4);

line([x_min*2, x_max*2], [0, 0], 'Color', 'k', 'LineWidth', 0.8); hold on;
plot(x_fit, y_fit, 'k-', 'LineWidth', 1.2);
errorbar(bx, by, bstd, 'o', 'Color', 'k', 'MarkerSize', 4, ...
    'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'LineWidth', 1.0, 'CapSize', 3);

xlim([x_min*1.5, x_max*1.5]);
ylim([-.5,.5]);
yticks(-0.4:0.2:0.4);
xlabel('$\nabla^{2} T$ (10$^{-10}$ $^\circ\mathrm{C}$ m$^{-2}$)', 'Interpreter', 'latex');
yl = ylabel('WSC (N m$^{-3} \times 10^7$)', 'Interpreter', 'latex');
ax = gca;
set(ax, 'FontSize', 13, 'Box', 'on', 'TickDir', 'in', 'XMinorTick', 'on', 'YMinorTick', 'on');
text(0.05, 0.88, sprintf('s = %.2f', slope), 'Units', 'normalized', 'FontSize', 13, 'FontName', 'Helvetica');
title('f)', 'Units', 'normalized', 'Position', [0, 1], 'HorizontalAlignment', 'left', 'FontSize', 15);

% -------------------------------------------------------------------------
% Subplot 7: Atlantic WSC & Crosswind SST Grad Contour (Panel d)
% -------------------------------------------------------------------------
ss(7) = subplot('Position', [px(4) py(1) axw axh]);
contourf(lona, lata, wscmeana_highpassed, -5e-7:1e-9:5e-7, 'LineStyle', 'None');
hold on;
contour(lona, lata, cwsgamean_highpassed, [-2e-7, 0], 'LineColor', 'black', 'LineStyle', '--');
contour(lona, lata, cwsgamean_highpassed, 0:100:100, 'LineColor', 'black', 'LineWidth', 1.2);
contour(lona, lata, cwsgamean_highpassed, 0:2e-7:2e-7, 'LineColor', 'black');
cm = redblue(101);
colorbar();
colormap(gca, cm);
worldmap3Atl(2);
xlim([-60 10]);
ylim([-30 30]);
caxis([-3e-8, 3e-8]);
xticks([-60 -30 0]);
xticklabels({'60°„W','30°„W','0°„'});
yticks([flip(-1*[0 10 20 30]) [10 20 30]]);
yticklabels({'30°„S','20°„S','10°„S','0°„','10°„N','20°„N','30°„N'});
title('d)', 'Units', 'normalized', 'Position', [0, 1], 'HorizontalAlignment', 'left');
ax = gca;
set(ax, 'FontSize', 15, 'Color', 0.6*[1 1 1], 'TickDir', 'out');
ax.YAxis.TickLabelGapOffset = -2.5;
ax.XAxis.TickLabelGapOffset = -4;

% Get actual rendered position of top panel 'd'
drawnow;
pos7 = ss(7).Position;

% -------------------------------------------------------------------------
% Subplot 8: Atlantic Crosswind SST Grad vs Curl Binned Plot (Panel h)
% -------------------------------------------------------------------------
ss(8) = subplot('Position', [pos7(1) py(2) pos7(3) axh1]);
x_min = -0.1; x_max = 0.1;
[bx, by, bstd, slope, x_fit, y_fit] = compute_binned_stats( ...
    cwsgamean_highpassed .* 1e5, wscmeana_highpassed .* 1e7, x_min, x_max, 20, 4);

line([x_min*2, x_max*2], [0, 0], 'Color', 'k', 'LineWidth', 0.8); hold on;
plot(x_fit, y_fit, 'k-', 'LineWidth', 1.2);
errorbar(bx, by, bstd, 'o', 'Color', 'k', 'MarkerSize', 4, ...
    'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'LineWidth', 1.0, 'CapSize', 3);

xlim([x_min*1.5, x_max*1.5]);
ylim([-.35,.35]);
yticks(-0.3:0.1:0.3);
xlabel('crosswind $\nabla T$ ($10^{-5}$ $^\circ\mathrm{C}$ m$^{-1}$)', 'Interpreter', 'latex');
ylabel('WSC (N m$^{-3} \times 10^7$)', 'Interpreter', 'latex');
ax = gca;
set(ax, 'FontSize', 13, 'Box', 'on', 'TickDir', 'in', 'XMinorTick', 'on', 'YMinorTick', 'on');
text(0.05, 0.88, sprintf('s = %.2f', slope), 'Units', 'normalized', 'FontSize', 13, 'FontName', 'Helvetica');
title('h)', 'Units', 'normalized', 'Position', [0, 1], 'HorizontalAlignment', 'left', 'FontSize', 15);

saveas(gcf,'test1.png');
close;