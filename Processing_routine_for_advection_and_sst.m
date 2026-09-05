longitude_p=ncread('CMEMS_GeostrophicU_Pacific.nc','longitude');
latitude_p=ncread('CMEMS_GeostrophicU_Pacific.nc','latitude');
u_p=ncread('CMEMS_GeostrophicU_Pacific.nc','ugos');
time=ncread('CMEMS_GeostrophicU_Pacific.nc','time');

longitude_p=longitude_p+360;

[I,J,K]=size(u_p);
u_p1=movmean(u_p,4,1,'omitnan');
u_p=movmean(u_p1,4,2,'omitnan');
clear u_p1

u_p = permute(u_p,[2 1 3]);
% 
% u_mean = movmean(u_p,360,3);
% 
% u_lowpassed1=movmean(u_mean,16,1,'omitnan');
% u_lowpassed2=movmean(u_lowpassed1,16,2,'omitnan');clear u_lowpassed1
% u_lowpassed3=movmean(u_mean-u_lowpassed2,16,1,'omitnan');
% u_lowpassed4=movmean(u_lowpassed3,16,2,'omitnan');clear u_lowpassed3
% up_highpassed=u_mean-(u_lowpassed4+u_lowpassed2);
% 
% up_highpassed=smoother_2d(up_highpassed,4);


load('IfremerdailyTAUandSSTAllPacific.mat');
clear taux tauy

sst=sst-273.15;

sst1=movmean(sst,4,1,'omitnan');
sst=movmean(sst1,4,2,'omitnan');
clear sst1

lonp = lons;clear lons
latp = lats;clear lats

sst = permute(sst,[2,1,3]);
[sstx,ssty] = zh_grad2(sst,lonp,latp);
clear ssty;


sstmean=movmean(sst,360,3,'omitnan');

sst_lowpassed1=movmean(sstmean,16,1,'omitnan');
sst_lowpassed2=movmean(sst_lowpassed1,16,2,'omitnan');clear sst_lowpassed1
sst_lowpassed3=movmean(sstmean-sst_lowpassed2,16,1,'omitnan');
sst_lowpassed4=movmean(sst_lowpassed3,16,2,'omitnan');clear sst_lowpassed3
sstp_highpassed=sstmean-(sst_lowpassed4+sst_lowpassed2);

sstmean=nanmean(sst,3);clear sst

sst_lowpassed1=movmean(sstmean,16,1,'omitnan');
sst_lowpassed2=movmean(sst_lowpassed1,16,2,'omitnan');clear sst_lowpassed1
sst_lowpassed3=movmean(sstmean-sst_lowpassed2,16,1,'omitnan');
sst_lowpassed4=movmean(sst_lowpassed3,16,2,'omitnan');clear sst_lowpassed3 
sstp_mean_highpassed=sstmean-(sst_lowpassed4+sst_lowpassed2); clear sst_lowpassed2 sst_lowpassed4

% 2. Create 2D target grid
[LON1, LAT1] = meshgrid(longitude_p, latitude_p);

ntime = size(sstx, 3);
matrix2_on_grid1 = zeros(length(latitude_p), length(longitude_p), ntime);
sstmean2p = zeros(length(latitude_p), length(longitude_p), ntime);

% 4. Interpolate slice by slice
for k = 1:ntime
    matrix2_on_grid1(:,:,k) = interp2(lonp, latp, sstx(:,:,k), LON1, LAT1, 'linear');
    sstmean2p(:,:,k) = interp2(lonp, latp, sstp_highpassed(:,:,k), LON1, LAT1, 'linear');
end
sstp_mean_highpassed = interp2(lonp, latp, sstp_mean_highpassed, LON1, LAT1, 'linear');

clear sstx
advection = -1*u_p(:,:,1:1461).*matrix2_on_grid1;clear u_p matrix2_on_grid1

advectionmean=movmean(advection,360,3,'omitnan');clear sst

advection_lowpassed1=movmean(advectionmean,16,1,'omitnan');
advection_lowpassed2=movmean(advection_lowpassed1,16,2,'omitnan');clear advection_lowpassed1
advection_lowpassed3=movmean(advectionmean-advection_lowpassed2,16,1,'omitnan');
advection_lowpassed4=movmean(advection_lowpassed3,16,2,'omitnan');clear advection_lowpassed3
advectionp_highpassed=advectionmean-(advection_lowpassed4+advection_lowpassed2);


advectionmean=nanmean(advection,3);

advection_lowpassed1=movmean(advectionmean,16,1,'omitnan');
advection_lowpassed2=movmean(advection_lowpassed1,16,2,'omitnan');clear advection_lowpassed1
advection_lowpassed3=movmean(advectionmean-advection_lowpassed2,16,1,'omitnan');
advection_lowpassed4=movmean(advection_lowpassed3,16,2,'omitnan');clear advection_lowpassed3
advectionmeanp_highpassed=advectionmean-(advection_lowpassed4+advection_lowpassed2);

[J,I]=size(advectionmeanp_highpassed);
A = find(longitude_p>260);
B = find(latitude_p>0);
for j=1:length(B)
    nanloc = find(advectionmeanp_highpassed(B(j),A)~=advectionmeanp_highpassed(B(j),A));
    advectionmeanp_highpassed(B(j),A(nanloc(1)):I)=0;
    sstp_mean_highpassed(B(j),A(nanloc(1)):I)=0;
end

load etopo_globe.mat
landmaskp =size(advectionmeanp_highpassed);
for i=1:I
    for j=1:J
        lonind = find(lon>=longitude_p(i));  
        latind = find(lat<=latitude_p(j));
        if topo(latind(1),lonind(1))<=0
            landmaskp(j,i)=topo(latind(1),lonind(1));
        else
            landmaskp(j,i)=NaN;
        end
    end
end    
clear topo;

for i =1:J
    nanlocy =  find(landmaskp(i,:)~=landmaskp(i,:));
    advectionmeanp_highpassed(i,nanlocy,:)=NaN;
    sstp_mean_highpassed(i,nanlocy,:)=NaN;
end

indsp = 1:120;
adv_regionp = squeeze(nanmean(advectionp_highpassed(:,indsp,:),2));
sst_regionp = squeeze(nanmean(sstmean2p(:,indsp,:),2));


%% ATLANTIC
longitude_a=ncread('CMEMS_GeostrophicU_Atlantic.nc','longitude');
latitude_a=ncread('CMEMS_GeostrophicU_Atlantic.nc','latitude');
u_a=ncread('CMEMS_GeostrophicU_Atlantic.nc','ugos');
time=ncread('CMEMS_GeostrophicU_Atlantic.nc','time');

[I,J,K]=size(u_a);
u_a1=movmean(u_a,4,1,'omitnan');
u_a=movmean(u_a1,4,2,'omitnan');
clear u_a1

u_a=permute(u_a,[2 1 3]);

% 
% u_mean = movmean(u_p,360,3);
% 
% u_lowpassed1=movmean(u_mean,16,1,'omitnan');
% u_lowpassed2=movmean(u_lowpassed1,16,2,'omitnan');clear u_lowpassed1
% u_lowpassed3=movmean(u_mean-u_lowpassed2,16,1,'omitnan');
% u_lowpassed4=movmean(u_lowpassed3,16,2,'omitnan');clear u_lowpassed3
% up_highpassed=u_mean-(u_lowpassed4+u_lowpassed2);
% 
% up_highpassed=smoother_2d(up_highpassed,4);


load('IfremerdailyWindstressAndSSTAtlantic.mat');
clear taux tauy

lona = lons;clear lons
lata = lats;clear lats

sst1=movmean(sst,4,1,'omitnan');
sst=movmean(sst1,4,2,'omitnan');
clear sst1

sst=sst-273.15;

sst = permute(sst,[2,1,3]);

[sstx,ssty] = zh_grad2(sst,lona+360,lata);
clear ssty;


sstmean=movmean(sst,360,3,'omitnan');

sst_lowpassed1=movmean(sstmean,16,1,'omitnan');
sst_lowpassed2=movmean(sst_lowpassed1,16,2,'omitnan');clear sst_lowpassed1
sst_lowpassed3=movmean(sstmean-sst_lowpassed2,16,1,'omitnan');
sst_lowpassed4=movmean(sst_lowpassed3,16,2,'omitnan');clear sst_lowpassed3
ssta_highpassed=sstmean-(sst_lowpassed4+sst_lowpassed2);

sstmean=nanmean(sst,3);clear sst

sst_lowpassed1=movmean(sstmean,16,1,'omitnan');
sst_lowpassed2=movmean(sst_lowpassed1,16,2,'omitnan');clear sst_lowpassed1
sst_lowpassed3=movmean(sstmean-sst_lowpassed2,16,1,'omitnan');
sst_lowpassed4=movmean(sst_lowpassed3,16,2,'omitnan');clear sst_lowpassed3 
ssta_mean_highpassed=sstmean-(sst_lowpassed4+sst_lowpassed2); clear sst_lowpassed2 sst_lowpassed4

% 2. Create 2D target grid
[LON1, LAT1] = meshgrid(longitude_a, latitude_a);

ntime = size(sstx, 3);
matrix2_on_grid1 = zeros(length(latitude_a), length(longitude_a), ntime);
sstmean2a = zeros(length(latitude_a), length(longitude_a), ntime);

% 4. Interpolate slice by slice
for k = 1:ntime
    matrix2_on_grid1(:,:,k) = interp2(lona, lata, sstx(:,:,k), LON1, LAT1, 'linear');
    sstmean2a(:,:,k) = interp2(lona, lata, ssta_highpassed(:,:,k), LON1, LAT1, 'linear');
end
ssta_mean_highpassed = interp2(lona, lata, ssta_mean_highpassed, LON1, LAT1, 'linear');

clear sstx
advection = -1*u_a(:,:,1:1461).*matrix2_on_grid1;clear u_p matrix2_on_grid1

advectionmean=movmean(advection,360,3,'omitnan');clear sst

advection_lowpassed1=movmean(advectionmean,16,1,'omitnan');
advection_lowpassed2=movmean(advection_lowpassed1,16,2,'omitnan');clear advection_lowpassed1
advection_lowpassed3=movmean(advectionmean-advection_lowpassed2,16,1,'omitnan');
advection_lowpassed4=movmean(advection_lowpassed3,16,2,'omitnan');clear advection_lowpassed3
advectiona_highpassed=advectionmean-(advection_lowpassed4+advection_lowpassed2);


advectionmean=nanmean(advection,3);

advection_lowpassed1=movmean(advectionmean,16,1,'omitnan');
advection_lowpassed2=movmean(advection_lowpassed1,16,2,'omitnan');clear advection_lowpassed1
advection_lowpassed3=movmean(advectionmean-advection_lowpassed2,16,1,'omitnan');
advection_lowpassed4=movmean(advection_lowpassed3,16,2,'omitnan');clear advection_lowpassed3
advectionmeana_highpassed=advectionmean-(advection_lowpassed4+advection_lowpassed2);


indsa = 161:321;
adv_regiona = squeeze(nanmean(advectiona_highpassed(:,indsa,:),2));
sst_regiona = squeeze(nanmean(sstmean2a(:,indsa,:),2));