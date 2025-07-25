# zamboni_et_al
Code for the analysis of 'Decoding injury responsive enhancers in the CNS for cell state targeting'

1. Pre-processing (getting data from raw sequencing files to Seurat objects, integrated and cleared of low quality cells and doublets, SCENICplus preprocessing, and CUTandTag) - Figure 1, Figure S1, Figure S11
   1. mgi_demultiplex.sh
   2. cellbender.sh
   3. archr_peak_calling.R
   4. create_object.R
   5. filtering.R
   6. merge_sample_objects.R
   7. config.yaml
   8. pycistopic_multiome.ipynb
   9. multiome_scenicplus.ipynb
   10. cutandtag_preprocess.R
      
2. Machine learning (prepare input data based on differentially accessible regions, train the cell type specific models and intepret their predictions) - Figure 2, 5, Figure S4, S5, S10, S11
   1. prepare_ml_input.R
   2. run_chrombpnet.sh
   3. peak_interpretation.ipynb
  
3. Explore (explorative analysis, such as differential expression and accessibility across cell types and injury timepoints, motif analysis, and processing and analysis of the enhancer aav reporter assay) - Figure 3, 4, 5, Figure S2, S3, S6, S7, S8, S9, S11, S12
   1. differential_analysis.R
   2. ml_interpretation.R
   3. explore_injury_differential_analysis.R
   4. ml_injury_interpretation.R
   5. explore_scenicplus.R
   6. enhancer_aav_reporter_assay.R
   7. cutandtag.R

