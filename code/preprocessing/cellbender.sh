# cellbender for ambient RNA removal

#pip install cellbender

directories=("ALL_1dpi_1" "ALL_1dpi_2" "ALL_1dpi_3" "ALL_28dpi_2" "ALL_28dpi_3"
             "ALL_3dpi_6" "ALL_7dpi_1" "ALL_7dpi_2" "ALL_7dpi_3" "ALL_U_1"
             "ALL_U_2" "ALL_U_3" "DEPL_1dpi_5" "DEPL_28dpi_5" "DEPL_3dpi_6"
             "DEPL_7dpi_5" "DEPL_U_5" "FT_1dpi_4" "FT_28dpi_4" "FT_7dpi_4"
             "FT_U_4")

for dir in "${directories[@]}"; do

cellbender remove-background \
--cuda \
--input $dir/outs/gex_raw_feature_bc_matrix.h5 \
--output $dir/outs/cellbended_gex_matrix_fpr001_seurat.h5 \
--fpr 0.001 \
--epochs 100;

done
