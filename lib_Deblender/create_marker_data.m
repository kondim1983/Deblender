function marker_data = create_marker_data(mixed_genes, mixed_data, marker_genes, marker_cell_type_index)

% % % mixed_genes: <column vector> all gene names (class: cell or double) of mixed dataset with 
% % % each each name corresponding to the row of "mixed_data"
% % % mixed_data: expression data matrix (class double) of mixed samples(non-logged)
% % % marker_genes: <column vector> genes (class: cell or double) with same nomenclature like "mixed_genes" 
% % % marker_cell_type_index: <column vector> (same length as marker_genes) with integer values 1...N where
% % % N is the number of cell types, which designates the source where each marker gene belongs to

% create marker data
Nc=max(marker_cell_type_index);
marker_data=[];
for i=1:1:Nc
    tf=ismember(marker_cell_type_index,i);
    source_markers=marker_genes(tf==1);
    tf1=ismember(mixed_genes, source_markers);
    marker_data_mean=mean(mixed_data(tf1==1,:),1);
    marker_data=[marker_data;marker_data_mean];
end

    
    






