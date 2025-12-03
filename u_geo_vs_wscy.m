longitude_p=ncread('CMEMS_GeostrophicU_Pacific.nc','longitude');
latitude_p=ncread('CMEMS_GeostrophicU_Pacific.nc','latitude');
u_p=ncread('CMEMS_GeostrophicU_Pacific.nc','ugos');
time=ncread('CMEMS_GeostrophicU_Pacific.nc','time');

u_lowpassed1=movmean(u_p,16,1,'omitnan');
u_lowpassed2=movmean(u_lowpassed1,16,2,'omitnan');clear u_lowpassed1
u_lowpassed3=movmean(u_p-u_lowpassed2,16,1,'omitnan');
u_lowpassed4=movmean(u_lowpassed3,16,2,'omitnan');clear u_lowpassed3
up_highpassed=u_p-(u_lowpassed4+u_lowpassed2);clear u_lowpassed2 u_lowpassed4

up_region=squeeze(nanmean(up_highpassed(:,1:121,:),2));
up_region_lp=movmean(up_region,360,2,'omitnan');


contourf(1:1462,latitude_p,up_region_lp,-1e-1:1e-3:1e-1,'LineStyle','None','edgecolor','None');hold on;
cm = redblue(101);
colorbar();
colormap(cm);
ylim([-30 0]);
caxis([-1e-1,1e-1])
contour(1:1462,latitude_p,up_region_lp,[0,1e-2,2e-2],'LineColor','black');
title(strcat('u_{geostrophic} 2005-2009'));

plot(nanmean(up_region_lp,2),latitude_p);

save('geostrophic_velocity_spatially_highpassed.mat','up_highpassed','longitude_p','latitude_p','time');


longitude_a=ncread('CMEMS_GeostrophicU_Atlantic.nc','longitude');
latitude_a=ncread('CMEMS_GeostrophicU_Atlantic.nc','latitude');
u_a=ncread('CMEMS_GeostrophicU_Atlantic.nc','ugos');
time=ncread('CMEMS_GeostrophicU_Atlantic.nc','time');

u_a=permute(u_a,[2 1 3]);

u_lowpassed1=movmean(u_a,16,1,'omitnan');
u_lowpassed2=movmean(u_lowpassed1,16,2,'omitnan');clear u_lowpassed1
u_lowpassed3=movmean(u_a-u_lowpassed2,16,1,'omitnan');
u_lowpassed4=movmean(u_lowpassed3,16,2,'omitnan');clear u_lowpassed3
ua_highpassed=u_a-(u_lowpassed4+u_lowpassed2);

up_region=squeeze(nanmean(ua_highpassed(:,161:321,:),2));
up_region_lp=movmean(up_region,360,2,'omitnan');


contourf(1:1462,latitude_a,up_region_lp,-1e-1:1e-3:1e-1,'LineStyle','None','edgecolor','None');hold on;
cm = redblue(101);
colorbar();
colormap(cm);
ylim([-30 30]);
caxis([-1e-1,1e-1])
contour(1:1462,latitude_a,up_region_lp,[0,1e-2,2e-2],'LineColor','black');
title(strcat('u_{geostrophic} 2005-2009'));

save('geostrophic_velocity_spatially_highpassed_ATL.mat','ua_highpassed','longitude_a','latitude_a','time');





load('wscy_over_beta.mat');
date_wscy=datestr(datenum('1-1-1')+time);
wscyQuik_lowpassed1=movmean(wscyQuik,16,1,'omitnan');
wscyQuik_lowpassed2=movmean(wscyQuik_lowpassed1,16,2,'omitnan');clear wscyQuik_lowpassed1
wscyQuik_lowpassed3=movmean(wscyQuik-wscyQuik_lowpassed2,16,1,'omitnan');
wscyQuik_lowpassed4=movmean(wscyQuik_lowpassed3,16,2,'omitnan');clear wscyQuik_lowpassed3
wscyQuik_highpassed=wscyQuik-(wscyQuik_lowpassed4+wscyQuik_lowpassed2);clear wscyQuik_lowpassed2 wscyQuik_lowpassed4
wscyQuik_region=squeeze(nanmean(wscyQuik_highpassed(:,361:521,:),2));
wscyQuik_region_lp=movmean(wscyQuik_region,360,2,'omitnan');

load('geostrophic_velocity_spatially_highpassed.mat');
date_ug=datestr(datenum('1970-1-1')+time/(24*3600));
up_region=squeeze(nanmean(up_highpassed(:,1:121,:),2));
up_region_lp=movmean(up_region,360,2,'omitnan');

load('wscy_over_beta_ATL.mat');
date_wscy_ATL=datestr(datenum('1-1-1')+time);
wscyQuik_lowpassed1_ATL=movmean(wscyQuik_A,16,1,'omitnan');
wscyQuik_lowpassed2_ATL=movmean(wscyQuik_lowpassed1_ATL,16,2,'omitnan');clear wscyQuik_lowpassed1_ATL
wscyQuik_lowpassed3_ATL=movmean(wscyQuik_A-wscyQuik_lowpassed2_ATL,16,1,'omitnan');
wscyQuik_lowpassed4_ATL=movmean(wscyQuik_lowpassed3_ATL,16,2,'omitnan');clear wscyQuik_lowpassed3_ATL
wscyQuik_highpassed_ATL=wscyQuik_A-(wscyQuik_lowpassed4_ATL+wscyQuik_lowpassed2_ATL);clear wscyQuik_lowpassed2_ATL wscyQuik_lowpassed4_ATL
wscyQuik_region_ATL=squeeze(nanmean(wscyQuik_highpassed_ATL(:,161:321,:),2));clear wscyQuik_A
wscyQuik_region_lp_ATL=movmean(wscyQuik_region_ATL,360,2,'omitnan');

load('geostrophic_velocity_spatially_highpassed_ATL.mat');
date_ug_ATL=datestr(datenum('1970-1-1')+time/(24*3600));
up_region_ATL=squeeze(nanmean(ua_highpassed(:,161:321,:),2));
up_region_lp_ATL=movmean(up_region_ATL,360,2,'omitnan');




plot(nanmean(wscyQuik_region_lp,2),latp);
hold on;
plot(nanmean(up_region_lp,2),latitude_p);
legend('$\frac{1}{\beta} \frac{\partial \left(\mathrm{\nabla}\times\vec{\tau}\right)\  }{\partial y}$','$u_{geostrophic}$','Interpreter','latex','fontsize',15);
grid on;
xlim([-0.04,0.04]);


%%
marg_h =[0.06 0.05]; 
marg_w =[0.05 0.05];
Nh = 2;
Nw = 5;
gap=[0.1 0.04];

axh = (1-sum(marg_h)-(Nh-1)*gap(1))/Nh; 
axw = (1-sum(marg_w)-2*gap(2))/2.5;

axh1 = (1-sum(marg_h)-(Nh-1)*gap(1))/Nh; 
axw1 = (1-sum(marg_w)-2*gap(2))/5;

px = [marg_w(1) marg_w(1)+axw+gap(2) marg_w(1)+axw*2+2*gap(2)];
py = [1-marg_h(1)-axh 1-marg_h(1)-axh*2-gap(1)]; 
set(gcf,'units','centimeters','position',[0,0,38.5,21]);
%%Pacific plot
ss(1)=subplot('Position',[px(1) py(1) axw axh]);
ax1=gca;
contourf(1:1461,latp,movmean(1e5*wscyQuik_region_lp/1025,4,1,'omitnan'),-2:1e-1:2,'LineStyle','None','edgecolor','None');hold on;
cm = redblue(101);
colorbar();
colormap(cm);
xticks([44,44+365,44+365*2,45+365*3]);
xticklabels({'2006','2007','2008','2009'});
ax1.XAxis.TickLabelGapOffset = -6;
ylim([-30 30]);
yticks([flip(-1*[0 10 20 30]) [10 20 30]]);
yticklabels({'30°„S','20°„S','10°„S','0°„','10°„N','20°„N','30°„N'});
caxis([1e5*-2e-5,1e5*2e-5])
contour(1:1461,latp,movmean(1e5*wscyQuik_region_lp/1025,4,1,'omitnan'),[0,1e5*5e-6,1e5*10e-6,1e5*10e-5],'LineColor','black');
ax1.XMinorTick='on';
ax1.XAxis.MinorTickValues=linspace(13,3373,113);
ttl=title('a) (${1}/{\beta \rho}) \partial [\mathrm{\nabla}\times\vec{\tau}] / \partial y$      ($\times$ 10$^{-5}$ ms$^{-1}$)','Interpreter','latex','Units','normalized','Position',[0.4,1.03]);
set(gca,'FontSize', 17,'TickDir','out'); 

ss(2)=subplot('Position',[px(2) py(1) axw axh]);
ax1=gca;
contourf(1:1462,latitude_p,1e2*up_region_lp,1e2*-2e-1:1e2*5e-3:1e2*2e-1,'LineStyle','None','edgecolor','None');hold on;
cm = redblue(101);
colorbar();
colormap(cm);
xticks([46,46+365,46+365*2,47+365*3]);
xticklabels({'2006','2007','2008','2009'});
ax1.XAxis.TickLabelGapOffset = -6;
ylim([-30 30]);
yticks([flip(-1*[0 10 20 30]) [10 20 30]]);
yticklabels({'30°„S','20°„S','10°„S','0°„','10°„N','20°„N','30°„N'});
caxis([-1e2*1e-1,1e2*1e-1])
contour(1:1462,latitude_p,up_region_lp,[0,1e2*5e-3,1e2*10e-2],'LineColor','black');
ax1.XMinorTick='on';
ax1.XAxis.MinorTickValues=linspace(15,3375,113);
ttl=title('b) Geostrophic velocity (cm/s)','Units','normalized','Position',[0.3,1.03]);
set(gca,'FontSize', 17,'TickDir','out'); 

ss(3)=subplot('Position',[px(3) py(1) axw1 axh1]);
ax1 = gca;
h1 = plot(ax1,movmean(1e6*nanmean(wscyQuik_region_lp,2)/1025,4,1,'omitnan'),latp, 'LineWidth', 1, 'color', 'blue');
hold on;
% Dummy plot for legend entry (invisible on ax1)
h3 = plot(ax1,nan, nan, 'LineWidth', 1, 'color', 'red');
grid on;
xlim([1e6*-0.000015,1e6*0.000015])
xlabel('$({1}/{\beta \rho}) \partial [\mathrm{\nabla}\times\vec{\tau}] / \partial y $     ($\times$ 10$^{-6}$ ms$^{-1}$)','Interpreter','latex','FontSize',15);
yticks([flip(-1*[0 10 20 30]) [10 20 30]]);
yticklabels({'30°„S','20°„S','10°„S','0°„','10°„N','20°„N','30°„N'});
set(gca,'FontSize', 16,'TickDir','out'); 
ax1.XLabel.Position(2) = ax1.XLabel.Position(2) +3.5; 
ax1.XLabel.Position(1) = ax1.XLabel.Position(1); 
ax1.XAxis.TickLabelGapOffset = -6;  % For bottom x-axis (blue)
ax1.XColor='blue';
% Secondary axes for wscy_smeridian (top x-axis)
ax2 = axes('Position', ax1.Position, ...
           'XAxisLocation', 'top', ...
           'YAxisLocation', 'left', ...
           'Color', 'none', ...
           'YColor', 'none', ...
           'XColor', 'red'); % Match u_g line color
hold(ax2,'on');
plot(ax2, 100*nanmean(up_region_lp,2),latitude_p, 'LineWidth', 1, 'color', 'red');
xlim(ax2, [-8,8]);
xlabel(ax2,'Geostrophic velocity (cm/s)')
ax2.XLabel.Position(2) = ax2.XLabel.Position(2)-0.8; 
ax2.XLabel.Position(1) = ax2.XLabel.Position(1)+1.8; 
ax2.XAxis.TickLabelGapOffset = -5;  % For top x-axis (orange)
%legend([h1, h3], {'$\frac{1}{\beta} \frac{\partial \left(\mathrm{\nabla}\times\vec{\tau}\right)\  }{\partial y}$','$u_{geostrophic}$'},'Interpreter','latex','fontsize',15);
set(gca,'FontSize', 16,'TickDir','out'); 
title('c)','Units','normalized','Position',[0,1.03],'FontSize',17);

%%ATLANTIC PLOT

ss(4)=subplot('Position',[px(1) py(2) axw axh]);
ax1=gca;
contourf(1:1461,lata,movmean(1e5*wscyQuik_region_lp_ATL/1025,4,1,'omitnan'),-2:1e-1:2,'LineStyle','None','edgecolor','None');hold on;
cm = redblue(101);
colorbar();
colormap(cm);
xticks([44,44+365,44+365*2,45+365*3]);
xticklabels({'2006','2007','2008','2009'});
ax1.XAxis.TickLabelGapOffset = -6;
ylim([-30 30]);
yticks([flip(-1*[0 10 20 30]) [10 20 30]]);
yticklabels({'30°„S','20°„S','10°„S','0°„','10°„N','20°„N','30°„N'});
caxis([1e5*-2e-5,1e5*2e-5])
contour(1:1461,lata,movmean(1e5*wscyQuik_region_lp_ATL/1025,4,1,'omitnan'),[0,1e5*5e-6,1e5*10e-6,1e5*10e-5],'LineColor','black');
ax1.XMinorTick='on';
ax1.XAxis.MinorTickValues=linspace(13,3373,113);
ttl=title('d) (${1}/{\beta \rho}) \partial [\mathrm{\nabla}\times\vec{\tau}] / \partial y$      ($\times$ 10$^{-5}$ ms$^{-1}$)','Interpreter','latex','Units','normalized','Position',[0.4,1.03]);
set(gca,'FontSize', 17,'TickDir','out'); 

ss(5)=subplot('Position',[px(2) py(2) axw axh]);
ax1=gca;
contourf(1:1462,latitude_a,1e2*up_region_lp_ATL,1e2*-2e-1:1e2*5e-3:1e2*2e-1,'LineStyle','None','edgecolor','None');hold on;
cm = redblue(101);
colorbar();
colormap(cm);
xticks([46,46+365,46+365*2,47+365*3]);
xticklabels({'2006','2007','2008','2009'});
ax1.XAxis.TickLabelGapOffset = -6;
ylim([-30 30]);
yticks([flip(-1*[0 10 20 30]) [10 20 30]]);
yticklabels({'30°„S','20°„S','10°„S','0°„','10°„N','20°„N','30°„N'});
caxis([-1e2*1e-1,1e2*1e-1])
contour(1:1462,latitude_a,up_region_lp_ATL,[0,1e2*5e-3,1e2*10e-2],'LineColor','black');
ax1.XMinorTick='on';
ax1.XAxis.MinorTickValues=linspace(15,3375,113);
ttl=title('e) Geostrophic velocity (cm/s)','Units','normalized','Position',[0.3,1.03]);
set(gca,'FontSize', 17,'TickDir','out'); 

ss(6)=subplot('Position',[px(3) py(2) axw1 axh1]);
ax1 = gca;
h1 = plot(ax1,movmean(1e6*nanmean(wscyQuik_region_lp_ATL,2)/1025,4,1,'omitnan'),lata, 'LineWidth', 1, 'color', 'blue');
hold on;
% Dummy plot for legend entry (invisible on ax1)
h3 = plot(ax1,nan, nan, 'LineWidth', 1, 'color', 'red');
grid on;
xlim([1e6*-0.000010,1e6*0.000010]);
xlabel('$({1}/{\beta \rho}) \partial [\mathrm{\nabla}\times\vec{\tau}] / \partial y $     ($\times$ 10$^{-6}$ ms$^{-1}$)','Interpreter','latex','FontSize',15);
yticks([flip(-1*[0 10 20 30]) [10 20 30]]);
yticklabels({'30°„S','20°„S','10°„S','0°„','10°„N','20°„N','30°„N'});
set(gca,'FontSize', 16,'TickDir','out'); 
ax1.XLabel.Position(2) = ax1.XLabel.Position(2)+4.1; 
ax1.XLabel.Position(1) = ax1.XLabel.Position(1); 
ax1.XAxis.TickLabelGapOffset = -7;  % For bottom x-axis (blue)
ax1.XColor='blue';
% Secondary axes for wscy_smeridian (top x-axis)
ax2 = axes('Position', ax1.Position, ...
           'XAxisLocation', 'top', ...
           'YAxisLocation', 'left', ...
           'Color', 'none', ...
           'YColor', 'none', ...
           'XColor', 'red'); % Match u_g line color
hold(ax2,'on');
plot(ax2, 100*nanmean(up_region_lp_ATL,2),latitude_a, 'LineWidth', 1, 'color', 'red');
xlim(ax2, [-6,6]);
xlabel(ax2,'Geostrophic velocity (cm/s)');
ax2.XLabel.Position(2) = ax2.XLabel.Position(2)-1; 
ax2.XLabel.Position(1) = ax2.XLabel.Position(1)+1.8; 
ax2.XAxis.TickLabelGapOffset = -5;  % For top x-axis (orange)
%legend([h1, h3], {'$\frac{1}{\beta} \frac{\partial \left(\mathrm{\nabla}\times\vec{\tau}\right)\  }{\partial y}$','$u_{geostrophic}$'},'Interpreter','latex','fontsize',15);
set(gca,'FontSize', 16,'TickDir','out'); 
title('f)','Units','normalized','Position',[0,1.03],'FontSize',17);