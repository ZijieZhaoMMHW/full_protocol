function [coef_MLR_RAW,r_raw,coef_LR_DETREND,coef_MLR_DETREND,r_detrend,...
    res_detrend,hc_t,t_full]=combining_together(ssta_raw,ssta_detrend,enso_line,...
    amo_line,pdo_line,ew_index)
%% Protocol

coef_MLR_RAW=NaN(size(ssta_raw,1),size(ssta_raw,2),4);
r_raw=NaN(size(ssta_raw,1),size(ssta_raw,2));

coef_LR_DETREND=NaN(size(ssta_raw,1),size(ssta_raw,2),4);

coef_MLR_DETREND=NaN(size(ssta_raw,1),size(ssta_raw,2),4);
r_detrend=NaN(size(ssta_raw,1),size(ssta_raw,2));
res_detrend=NaN(size(ssta_raw,1),size(ssta_raw,2),size(ssta_detrend,3));

hc_t=NaN(size(ssta_detrend,1),size(ssta_detrend,2),5);
t_full=(1:size(ssta_detrend,3))';

parfor i=1:size(ssta_raw,1);
    tic
    coef1_here=NaN(size(ssta_raw,2),4);
    r1_here=NaN(size(ssta_raw,2),1);
    
    coef2_here=NaN(size(ssta_raw,2),4);
    
    coef3_here=NaN(size(ssta_raw,2),4);
    r3_here=NaN(size(ssta_raw,2),1);
    res_here=NaN(size(ssta_raw,2),size(ssta_detrend,3));
    
    hc_here=NaN(size(ssta_detrend,2),5);
    
    for j=1:size(ssta_raw,2);
        
        % Part 1
        
        fixed_here=squeeze(ssta_raw(i,j,:));
        
        if nansum(isnan(fixed_here))==0 && length(unique(fixed_here))~=1
           
           data_used=[ones(size(enso_line,1),1) (1:size(enso_line,1))' enso_line(:,end) amo_line(:,end) pdo_line(:,end)];
            
           
           [b,bint,r,~,stats] = regress(fixed_here,data_used);
           
           coef1_here(j,1)=b(2);
           coef1_here(j,2)=b(3);
           coef1_here(j,3)=b(4);
           coef1_here(j,4)=b(5);
           
           r1_here(j)=stats(1);
        end
        
        % Part 2
        
        fixed_here=squeeze(ew_index(i,j,:));
        
        if nansum(isnan(fixed_here))==0 && nansum(fixed_here==0)~=length(fixed_here) && length(unique(fixed_here))~=1
            
            data_used=[(1:size(enso_line,1))' enso_line(:,end) amo_line(:,end) pdo_line(:,end)];
            
           fixed_bin=(double(fixed_here~=0));
           fixed_bin_cell=cell(length(fixed_bin),1);
           fixed_bin_cell(fixed_bin==0)={'Normal'};
           fixed_bin_cell(fixed_bin==1)={'Hot'};
           fixed_bin_cell=categorical(fixed_bin_cell);
           
           
           [B_fixed,dev_fixed,stats_fixed] = mnrfit(data_used,fixed_bin_cell);
           
           coef2_here(j,1)=B_fixed(2);
           coef2_here(j,2)=B_fixed(3);
           coef2_here(j,3)=B_fixed(4);
           coef2_here(j,4)=B_fixed(5);
           
        end
        
        % Part 3
        
        fixed_here=squeeze(ssta_detrend(i,j,:));
        
        if nansum(isnan(fixed_here))==0 && length(unique(fixed_here))~=1
           
           data_used=[ones(size(enso_line,1),1) (1:size(enso_line,1))' enso_line(:,end) amo_line(:,end) pdo_line(:,end)];
            
           
           [b,bint,r,~,stats] = regress(fixed_here,data_used);
           
           coef3_here(j,1)=b(2);
           coef3_here(j,2)=b(3);
           coef3_here(j,3)=b(4);
           coef3_here(j,4)=b(5);
           
           r3_here(j)=stats(1);
           
           res_here(j,:)=r;
        end
        
        % Part 4
        
        ew_here=squeeze(ew_index(i,j,:));
        r_here=squeeze(res_here(j,:));
        
        ct_here=nanmean(squeeze(coef3_here(j,1)).*t_full(ew_here==1));
        censo_here=nanmean(squeeze(coef3_here(j,2)).*enso_line(ew_here==1,end));
        camo_here=nanmean(squeeze(coef3_here(j,3)).*amo_line(ew_here==1,end));
        cpdo_here=nanmean(squeeze(coef3_here(j,4)).*pdo_line(ew_here==1,end));
        
        hc_here(j,1)=nanmean(ct_here);
        hc_here(j,2)=nanmean(censo_here);
        hc_here(j,3)=nanmean(camo_here);
        hc_here(j,4)=nanmean(cpdo_here);
        
        hc_here(j,5)=nanmean(r_here(ew_here==1));
        
    end
    
    coef_MLR_RAW(i,:,:)=coef1_here;
    r_raw(i,:)=r1_here;
    
    coef_LR_DETREND(i,:,:)=coef2_here;
    
    coef_MLR_DETREND(i,:,:)=coef3_here;
    r_detrend(i,:)=r3_here;
    res_detrend(i,:,:)=res_here;
    
    hc_t(i,:,:)=hc_here;
    
    toc
end