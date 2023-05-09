clear all;clc;close all
%8bit files pre-processed gaussian blur and background subtracted
folder_path= uigetdir();
file_finder=dir([folder_path,'*_g_bc*tif*']);
%parameters of the simulation, estimated by running a few scenarios
thresh_im=5;
dist_merge=80;
noise_thresh=50;
pixel_2_mic=0.455;
%set criteria for identifying a seed
c_seeds=[];
nog3_seeds=[];
seed_min=0.25;
seed_max=0.75;
for ff=1:length(file_finder)
    struct_ff=file_finder(ff);
    file_to_read=struct_ff.name
    file_read=imread(strcat(folder_path,file_to_read));
    %notochord masks are created and stored 
    th_nm=strcat(folder_path,'bn_files/',extractBefore(file_to_read,'.tif'),'_chord_mask.txt');
    if(~isfile(th_nm))
        figure(1);imshow(imcomplement(imadjust(file_read)))
        hold on
        r2=drawpolygon;
        c3=r2.Position(:,1);
        c4=r2.Position(:,2);
        bw1=roipoly(file_read,r2.Position(:,1),r2.Position(:,2));
        imshow(bw1)
        writematrix(bw1,strcat(folder_path,'bn_files/',extractBefore(file_to_read,'.tif'),'_chord_mask.txt'));
        clf(1)
    else
        bw1=load(th_nm);
    end
    noto_region=uint8(bw1).*file_read;
    bn_noto=noto_region>thresh_im;
    %connect regions and label them using bwconncomp & labelmatrix
    CC=bwconncomp(bn_noto);
    labels=labelmatrix(CC);
    label_colors=label2rgb(labels,"jet");
    %use regionprops to estimate area and coordinate list
    props=regionprops(labels,{"Area","PixelIdxList","Centroid"});
    area_first=cell2mat({props(:).Area});
    area_indices=find(area_first<=noise_thresh);
    props_cor=props;
    props_cor(area_indices)=[];
    %props sort by Centroid Col
    all_centroids=cell2mat(reshape({props_cor(:).Centroid},size(props_cor,1),1));
    all_cols=all_centroids(:,1);
    [all_cols_sorted,indices]=sort(all_cols,'ascend');
    props_sorted=props_cor(indices);
    %merge similar groups together
    groups_to_merge=join_peaks_close(all_cols_sorted,dist_merge);
    %create new groups based on merging and also confirm centroid positions
    new_area_list=zeros(length(groups_to_merge),1);
    new_ap_extent=zeros(length(groups_to_merge),1);
    new_dv_extent=zeros(length(groups_to_merge),1);
    new_area_pos=zeros(length(groups_to_merge),2);
    for gg = 1:length(groups_to_merge)
        list_gg=groups_to_merge{gg};
        new_area_list(gg)=sum([props_sorted(list_gg).Area]);
        all_cell_px={props_sorted(list_gg).PixelIdxList};
        comb_all_cell_px=unique(cell2mat(cat(1,all_cell_px')));
        [r,c]=ind2sub(size(labels),comb_all_cell_px);
        m1=mean(r);
        m2=mean(c);
        new_area_pos(gg,1)=new_area_list(gg)*pixel_2_mic*pixel_2_mic;
        new_area_pos(gg,2)=m2*pixel_2_mic;
        new_ap_extent(gg,1)=(max(c)-min(c))*pixel_2_mic;
        %find DV extent
        bn_each_col=noto_region(:,round(m2));
        non_zero_all=find(bn_each_col);
        min_zero=min(non_zero_all);
        max_zero=max(non_zero_all);
        new_dv_extent(gg,1)=(max_zero-m1)/(max_zero-min_zero);
        if(contains(file_to_read,'Control'))
            c_seeds=[c_seeds;new_dv_extent(gg,1)];
        else
            nog3_seeds=[nog3_seeds;new_dv_extent(gg,1)];
        end        
    end
    if(contains(file_to_read,'Control'))
        figure(1)
        h1=scatter(new_area_pos(:,2),new_area_pos(:,1),50,'MarkerFaceColor','#77AADD','MarkerEdgeColor','#77AADD');
        hold on
        %
        figure(2)
        h1p=scatter(new_area_pos(:,2),new_ap_extent(:,1),50,'MarkerFaceColor','#77AADD','MarkerEdgeColor','#77AADD');
        hold on
        %
        figure(3)
        h1pp=scatter(new_area_pos(:,2),new_dv_extent(:,1),50,'MarkerFaceColor','#77AADD','MarkerEdgeColor','#77AADD');
        hold on
        
    else
        figure(1)
        h2=scatter(new_area_pos(:,2),new_area_pos(:,1),50,'MarkerFaceColor','#FFAABB','MarkerEdgeColor','#FFAABB');
        %
        figure(2)
        h2p=scatter(new_area_pos(:,2),new_ap_extent(:,1),50,'MarkerFaceColor','#FFAABB','MarkerEdgeColor','#FFAABB');
        %
        figure(3)
        h2pp=scatter(new_area_pos(:,2),new_dv_extent(:,1),50,'MarkerFaceColor','#FFAABB','MarkerEdgeColor','#FFAABB');
        
    end
end    
%properties
%figure1
f=figure(1);
f.Position=[100 100 950 600]
box on;
legend([h1,h2],'CONTROL','NOG3')
xlabel('AP Axis (\mum)','interpreter','tex')
ylabel('Area of TP1 Clusters (\mum^2)', 'interpreter', 'tex')
set(gca,'fontsize',24,'fontname','ariel','linewidth',2)
%figure 2
f=figure(2);
f.Position=[100 100 950 600]
xlim([0 2600])
box on;
legend([h1p,h2p],'CONTROL','NOG3')
xlabel('AP Axis (\mum)','interpreter','tex')
ylabel('Cluster Extent Along AP (\mum)','interpreter','tex')
set(gca,'fontsize',24,'fontname','ariel','linewidth',2)
%figure 3
f=figure(3);
f.Position=[100 100 950 600]
x1=0:2600;
y1=seed_min*ones(length(x1),1);
y2=seed_max*ones(length(x1),1);
plot(x1,y1,'k--','linewidth',2)
plot(x1,y2,'k--','linewidth',2)
xlim([0 2600])
box on;
legend([h1pp,h2pp],'CONTROL','NOG3')
xlabel('AP Axis (\mum)','interpreter','tex')
ylabel('Cluster Extent Along DV')
set(gca,'fontsize',24,'fontname','ariel','linewidth',2)


%figure as scattergory
c_seeds_f1=c_seeds>0.25;
c_seeds_f2=c_seeds<0.75;
c_seeds_f=c_seeds(c_seeds_f1&c_seeds_f2);
c_frac=1-(length(c_seeds_f)/length(c_seeds))

nog3_seeds_f1=nog3_seeds>0.25;
nog3_seeds_f2=nog3_seeds<0.75;
nog3_seeds_f=nog3_seeds(nog3_seeds_f1 & nog3_seeds_f2);
nog3_frac=1-(length(nog3_seeds_f)/length(nog3_seeds))

