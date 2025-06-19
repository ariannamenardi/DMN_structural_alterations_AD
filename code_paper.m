
%% code for multiple linear regression models
% run separate for A+ and A- individuals

for i=1:size(DMN,2)
clear tbl lme* F n 

n=DMN(i)

tbl= table(zscore(N),zscore(sulc_depth(:,n)),zscore(thickness(:,n)), zscore(gyrification(:,n)),...
    zscore(tau_burden(:,i)),group, zscore(str2double(demographics.Age)),...
    normalize(MMSE), zscore(education),'VariableNames',...
    {'cognition','depth','thickness','gyrification','tau',...
    'group','age', 'MMSE', 'education'}); 


%check multicollinearity
clear x vif_value
x=cat(2,tbl.depth,tbl.thickness,tbl.gyrification,tbl.tau,tbl.age,tbl.MMSE,tbl.education);
vif_value(i,:)=vif(x);

lme1 = fitlm(tbl,['cognition~MMSE*(depth+thickness+gyrification)+age+education*(depth+thickness+gyrification)+tau*(depth+thickness+gyrification)']); 

outlier = find(isoutlier(lme1.Residuals.Raw)); %check potential outliers

lme2 = fitlm(tbl,'cognition~MMSE*(depth+thickness+gyrification)+age+education*(depth+thickness+gyrification)+tau*(depth+thickness+gyrification)', 'Exclude',[outlier]) %after correcting based on VIF


p_values_model(i)=lme2.ModelFitVsNullModel.Pvalue;
p_values_thick_mmse(i)=table2array(lme2.Coefficients("thickness:MMSE","pValue"));
p_values_depth_mmse(i)=table2array(lme2.Coefficients("depth:MMSE","pValue"));
p_values_gyri_mmse(i)=table2array(lme2.Coefficients("gyrification:MMSE","pValue"));

p_values_thick_tau(i)=table2array(lme2.Coefficients("thickness:tau","pValue"));
p_values_depth_tau(i)=table2array(lme2.Coefficients("depth:tau","pValue"));
p_values_gyri_tau(i)=table2array(lme2.Coefficients("gyrification:tau","pValue"));

p_values_thick_edu(i)=table2array(lme2.Coefficients("thickness:education","pValue"));
p_values_depth_edu(i)=table2array(lme2.Coefficients("depth:education","pValue"));
p_values_gyri_edu(i)=table2array(lme2.Coefficients("gyrification:education","pValue"));


%slope of significant interactions
data=plotInteraction(lme2,'MMSE','thickness','predictions'); % do the same for all morphology measures
set(gcf, 'Visible','off') %to turn off the automatic opening of the image
x_c=data(1,1).XData; 
y_c=data(1,1).YData;
m_thick_MMSElow(i) = (y_c(2)-y_c(1))/(x_c(2)-x_c(1)) ; %slope
x_m=data(2,1).XData; 
y_m=data(2,1).YData;
m_thick_MMSEmed(i) = (y_m(2)-y_m(1))/(x_m(2)-x_m(1)) ;
x_e=data(3,1).XData;
y_e=data(3,1).YData;
m_thick_MMSEhigh(i) = (y_e(2)-y_e(1))/(x_e(2)-x_e(1)) ;


%slope of significant interactions
data=plotInteraction(lme2,'tau','thickness','predictions'); % do the same for all morphology measures
set(gcf, 'Visible','off') 
x_c=data(1,1).XData; 
y_c=data(1,1).YData;
m_thick_taulow(i) = (y_c(2)-y_c(1))/(x_c(2)-x_c(1)) ; %slope
x_m=data(2,1).XData;
y_m=data(2,1).YData;
m_thick_taumed(i) = (y_m(2)-y_m(1))/(x_m(2)-x_m(1)) ; 
x_e=data(3,1).XData; 
y_e=data(3,1).YData;
m_thick_tauhigh(i) = (y_e(2)-y_e(1))/(x_e(2)-x_e(1)) ; 


%slope of significant interactions
data=plotInteraction(lme2,'education','thickness','predictions');% do the same for all morphology measures
set(gcf, 'Visible','off') 
x_c=data(1,1).XData; 
y_c=data(1,1).YData;
m_thick_lowedu(i) = (y_c(2)-y_c(1))/(x_c(2)-x_c(1)) ; %slope
x_m=data(2,1).XData; 
y_m=data(2,1).YData;
m_thick_mededu(i) = (y_m(2)-y_m(1))/(x_m(2)-x_m(1)) ; 
x_e=data(3,1).XData; 
y_e=data(3,1).YData;
m_thick_highedu(i) = (y_e(2)-y_e(1))/(x_e(2)-x_e(1)) ; 



figure;
 plotInteraction(lme2,'MMSE','depth','predictions') %interaction plots
 ax=gca;
 ax.FontSize=18;

 figure;
 plotInteraction(lme2,'tau','depth','predictions')
 ax=gca;
 ax.FontSize=18;

 figure;
 plotInteraction(lme2,'education','depth','predictions')
 ax=gca;
 ax.FontSize=14;


clear data x_* y_*
end

[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p_values_model); % FDR correction

%% Pearson's correlation between baseline structural measures and tau burden
% at follow-up

for i=1:length(DMN)
    n= DMN(i);
    [r_temp p_temp]=corr(zscore(sulc_depth(:,n)),tau_burden_y2(:,i),'rows','pairwise','Type','Pearson');
    r_y2(i,:)=r_temp;
    p_y2(i,:)=p_temp;
    clear n *_temp
end

[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p_y2); % FDR correction

%%  vertices correspondance
% Desikian
[vertices_rh,label_rh,colortable_rh]=read_annotation('/dati/sw/freesurfer/subjects/fsaverage5/label/rh.aparc.annot');
for i=1:length(vertices_rh)
[li_rh loc_rh]=find(colortable_rh.table==label_rh(i));
roi_DK_rh(i)=li_rh;
clear li_rh loc_rh
end
roi_names_rh=colortable_rh.struct_names(roi_DK_rh);

[vertices_lh,label_lh,colortable_lh]=read_annotation('/dati/sw/freesurfer/subjects/fsaverage5/label/lh.aparc.annot');
for i=1:length(vertices_lh)
[li_lh loc_lh]=find(colortable_lh.table==label_lh(i));
roi_DK_lh(i)=li_lh;
clear li_lh loc_lh
end
roi_names_lh=colortable_rh.struct_names(roi_DK_lh);

vertices_DK=cat(1,roi_DK_lh',roi_DK_rh');

%Schaefer
vertices_Schaefer=readtable('ENIGMA/matlab/shared/parcellations/schaefer_400_fsa5.csv');
vertices_Schaefer.Properties.VariableNames={'vertices_Schaefer'};

%combine to see the correpondence
vertices=cat(2,table2array(vertices_Schaefer),vertices_DK);
vertices_DK=array2table(vertices_DK);
vertices=cat(2,vertices_Schaefer,vertices_DK, cell2table(cat(1,roi_names_lh, roi_names_rh)));
vertices.Properties.VariableNames(3)={'DK_ROI'};
vertices=sortrows(vertices);

%% tau burden

Tau_PET_6months % TAU data provided by ADNI
for i=1:length(DMN)

    n=DMN(i)
    idx=find(vertices.vertices_Schaefer==n);
    region_DK=vertices.DK_ROI(idx);
    region_correspondance{i}=region_DK;

    [C, ~, ic] = unique(region_DK);
    counts = accumarray(ic(:), 1);
    [mm l]=max(counts);
    kk=find(ismember(region_DK,C(l)));
    max_region=mean(Tau_PET_6months(:,kk),2); %DK roi that corrisponds to Schaefer region with maximum vertex overlap

    clear counts
    nCats = numel(C);
        for q = 1:nCats
            idx = ic == q;
            means(q,:) = mean(Tau_PET_6months(:,idx),2);
            counts(q) = sum(idx);
            clear idx
            avg(q,:)=means(q,:) .* counts(q);
        end 
        
        weighted_mean(:,i) = sum(avg,1,'omitnan') / sum(counts); %weighted average based on number of vertices


        clear n idx region_DK C ic counts mm l kk avg means 
end