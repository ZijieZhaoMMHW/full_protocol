function [coef_MLR_RAW,r_raw,coef_LR_DETREND,coef_MLR_DETREND,r_detrend,...
    res_detrend,hc_t]=protocol(sst,year_range,clim_add,lat,enso_line,amo_line,pdo_line,...
    land_index)

%protocol - Execute all experiments in Liu et al. (in prep.)
%  Syntax
%
%  [coef_MLR_RAW,r_raw,coef_LR_DETREND,coef_MLR_DETREND,r_detrend,...
%    res_detrend,hc_t]=protocol(sst,year_range,clim_add,lat,enso_line,amo_line,pdo_line,...
%    land_index)
%
%  Description
%
%  [coef_MLR_RAW,r_raw,coef_LR_DETREND,coef_MLR_DETREND,r_detrend,...
%    res_detrend,hc_t]=protocol(sst,year_range,clim_add,lat,enso_line,amo_line,pdo_line,...
%    land_index) takes sea surface temperature (sst; x-y-t), the range of
%    year corresponding to sst (year_range; a two - element vector e.g.
%    [1957 2017], climate base line (clim_add; x-y-t), time series of ENSO
%    index (enso_line; t-1), AMO index (amo_line; t-1), PDO index
%    (pdo_line; t-1), and the label of land area (land_index, where the
%    land is masked as nan.
%
%  Output Arguments
%   
%   coef_MLR_RAW - A 3D matrix (x-y-4) containing coefficients from MLR
%   based on raw SST anomalies, the 4 layers separately correspond to
%   coefficients of T, ENSO, AMO, and PDO.
%   
%   r_raw - A 2D matrix (x-y) containing r-square from MLR based on raw SST
%   anomalies.
%
%   coef_LR_DETREND - A 3D matrix (x-y-4) containing coefficients from LR
%   based on detrend SST anomalies, the 4 layers separately correspond to
%   coefficients of T, ENSO, AMO, and PDO.
%
%   coef_MLR_DETREND - The same as coef_MLR_RAW, but for detrend SST
%   anomalies.
%
%   r_detrend - The same as r_raw, but for detrend SST anomalies.
%
%   res_detrend - A 3D matrix (x-y-t), containing residuals from the MLR
%   based on detrend SST.
%
%   hc_t - A 3D matrix (x-y=5), containing the heat contributions pf T,
%   ENSO, AMO, PDO, and Residuals from the MLR based on detrend SST.

m_y_used=[sort(repmat((year_range(1):year_range(2))',12,1)) repmat((1:12)',...
    year_range(2)-year_range(1)+1,1)];

ssta_raw=NaN(size(sst,1),size(sst,2),size(sst,3));
threshold_clim=NaN(size(sst,1),size(sst,2),size(sst,3));
for i=1:12
    index_here=m_y_used(:,2)==i;
    ssta_raw(:,:,index_here)=sst(:,:,index_here)-nanmean(sst(:,:,index_here),3);
    threshold_clim(:,:,index_here)=repmat(quantile(sst(:,:,index_here),0.9,3),1,1,nansum(index_here));
end

area_used=repmat(cosd(double(lat')),size(ssta_raw,1),1);
ssta_raw(isnan(ssta_raw))=0;
[eof_fixed,pc_fixed,expvar_fixed] = eof(ssta_raw.*repmat(area_used,1,1,size(ssta_raw,3)),'mask',~isnan(land_index));
t_used=NaN(size(pc_fixed,1),1);
for i=1:size(pc_fixed,1)
    pc_here=pc_fixed(i,:)./nanstd(pc_fixed(i,:));
    t_here=polyfit((1:size(pc_here,2))',pc_here(:),1);
    t_used(i)=t_here(1);
end
v=expvar_fixed./nansum(expvar_fixed);
p_index=1:length(t_used);
index_max=find(abs(t_used)==nanmax(abs(t_used)));

eof_fixed=eof_fixed(:,:,[index_max  p_index(p_index~=index_max)]);
pc_fixed=pc_fixed([index_max  p_index(p_index~=index_max)],:);
v=v([index_max  p_index(p_index~=index_max)]);

[a1,reg_a1,~,~,~] = Trend_Transfer((pc_fixed(~isnan(sum(pc_fixed,2)),:))',eof_fixed(:,:,~isnan(sum(pc_fixed,2))),v(~isnan(sum(pc_fixed,2))));

a1=a1';

ssta_detrend=ssta_raw.*repmat(area_used,1,1,size(ssta_raw,3))...
    -reof(reg_a1,a1,1);

ssta_detrend=ssta_detrend./repmat(area_used,1,1,size(ssta_raw,3));

sst_re=clim_add+ssta_detrend;

ew_index=NaN(size(ssta_detrend));

ew_index(sst_re>threshold_clim)=1;
ew_index(sst_re<threshold_clim)=0;

[coef_MLR_RAW,r_raw,coef_LR_DETREND,coef_MLR_DETREND,r_detrend,...
    res_detrend,hc_t,~]=combining_together(ssta_raw,ssta_detrend,enso_line,...
    amo_line,pdo_line,ew_index);