function groups_to_merge_l=join_peaks_close(all_cols_sorted,dist_merge)
    if(~iscolumn(all_cols_sorted))
        all_cols_sorted=all_cols_sorted';
    end    
    num_clusters=1:length(all_cols_sorted);
    list_pairs=[all_cols_sorted(1:end-1) all_cols_sorted(2:end)];
    row_wise_dist=abs(list_pairs(:,2)-list_pairs(:,1));
    find_breaks=find(row_wise_dist>dist_merge);
    groups_to_merge_l=[];
    ini_pos=1;
    for lc_i = 1:length(find_breaks)
        lc=find_breaks(lc_i);
        groups_to_merge_l=[groups_to_merge_l;{num_clusters(ini_pos:lc)}];
        ini_pos=lc+1;
    end    
    groups_to_merge_l=[groups_to_merge_l;{num_clusters(ini_pos:end)}];  
end