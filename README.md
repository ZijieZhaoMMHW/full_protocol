Examine the Influence of Climate Change on Extreme Warming's Response to Climate Modes
==================================================================

This repo contains part of codes to support what presented in Liu et al. (in prep.). Here I present the full experiment protocol and some test data from observations (SST from HADISST).

Data Description
-------------

`sst.mat` contains the sea surface temperature from HADISST during 1957 to 2017 in monthly resolution (360 * 180 * 732), `lon_lat_land_index` contains the longitude (360), latitude (180), and land index masked by NaN (360 * 180), `clim_add.mat` contains the base - line climatology (360 * 180 * 732), `clim_index.mat` contains Nino34, AMO, PDO indexes during 1957 to 2017 (732). 

Execute the Experiment
-------------

It is as easy as:

```
[coef_MLR_RAW,r_raw,coef_LR_DETREND,coef_MLR_DETREND,r_detrend,...
    res_detrend,hc_t]=protocol(sst_used,[1957 2017],clim_clim,lat,enso_line,amo_line,pdo_line,...
    land_index);
```

It returns several outputs, which are separately:

```
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
%   hc_t - A 3D matrix (x-y-5), containing the heat contributions for T,
%   ENSO, AMO, PDO, and Residuals from the MLR based on detrend SST.
```

MLR (Multiple Linear Regression) based on raw SST anomalies
-------------

![Image text](https://github.com/ZijieZhaoMMHW/full_protocol/blob/master/MLR_RAW.png)

MLR based on detrend SST anomalies
-------------

![Image text](https://github.com/ZijieZhaoMMHW/full_protocol/blob/master/MLR_DETREND.png)

LR (Logistic Regression) based on detrend SST anomalies
-------------

![Image text](https://github.com/ZijieZhaoMMHW/full_protocol/blob/master/LR_DETREND.png)

Heat Contribution from each covariate
-------------

![Image text](https://github.com/ZijieZhaoMMHW/full_protocol/blob/master/hc_try.png)
