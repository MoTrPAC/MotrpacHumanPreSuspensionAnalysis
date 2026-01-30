
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MotrpacHumanPreSuspensionAnalysis

<!-- badges: start -->
<!-- badges: end -->

# Installation

``` r
devtools::install_github("MoTrPAC/MotrpacHumanPreSuspensionAnalysis")
```

## Package Use and Structure

This package is the public release version of the results from the
Molecular Transducers of Physical Activity Consortium. The primary aim
of the package is to facilitate broad access to consortium-level
findings by distributing summary statistics across all reported
analyses, enabling reproducible downstream interpretation.

To protect participant privacy and comply with data-use governance
policies, individual-level (subject-level) molecular data are not
included in this package and are available only through formal data
access requests to the MoTrPAC consortium. Clinical endpoints, including
select clinical analyte panels, are available for all subjects.

### Clinical and phenotypic data

Clinical and phenotypic should be reference directly as data objects.
The `pheno` object is the main source for all of the phenotypic data,
and contains most experimental design parameters.

The package has data objects that come as a list with 2 items:

- `data` the raw data
- `dict` the dictionary for the data

``` r

library(MotrpacHumanPreSuspensionData)

pheno_data = load_pheno(load_acute_only = TRUE)
```

You can also directly call the `pheno` data object instead of using
`load_pheno` but beware that this will also load in the participants
that completed their chronic training, which are excluded from the vast
majority of our analysis.

Don’t waste time guessing what 1’s and 2’s mean. Attach your dictionary
to your data:

``` r

pheno_formatted <- attach_dictionary(pheno_data)

pheno_formatted[1:5, c("study", "sex_psca")]
#>                       study sex_psca
#> 11001010209 Adult Sedentary   Female
#> 11001010210 Adult Sedentary   Female
#> 11001010212 Adult Sedentary   Female
#> 11001010213 Adult Sedentary   Female
#> 11001010214 Adult Sedentary   Female
```

There are many different data items, if you directly look through each
of the data files, this can get very confusing.

``` r
grep("^cln", data(package = "MotrpacHumanPreSuspensionData")$results[, "Item"],
     value = TRUE)
#>  [1] "cln_chemistry_t02_plasma_lab_ck_results"               
#>  [2] "cln_chemistry_t02_plasma_lab_conv_results"             
#>  [3] "cln_chemistry_t02_plasma_lab_crt_results"              
#>  [4] "cln_chemistry_t02_plasma_lab_glc_results"              
#>  [5] "cln_chemistry_t02_plasma_lab_ins_results"              
#>  [6] "cln_curated_accel_derived_variables_baseline"          
#>  [7] "cln_curated_acute_bout"                                
#>  [8] "cln_curated_biospec"                                   
#>  [9] "cln_curated_screening"                                 
#> [10] "cln_raw_1rm_assess_test"                               
#> [11] "cln_raw_24_hr_food_record"                             
#> [12] "cln_raw_activity_monitor_record"                       
#> [13] "cln_raw_acute_endurance_ex_test"                       
#> [14] "cln_raw_acute_resistance_ex_test"                      
#> [15] "cln_raw_adi_collect_tp1_pre_ex"                        
#> [16] "cln_raw_adi_collect_tp2_24_hr_post_ex"                 
#> [17] "cln_raw_adi_collect_tp2_45_min_post_ex"                
#> [18] "cln_raw_adi_collect_tp2_4_hr_post_ex"                  
#> [19] "cln_raw_adverse_event_log"                             
#> [20] "cln_raw_analyticdatasets_biospecview"                  
#> [21] "cln_raw_analyticdatasets_key"                          
#> [22] "cln_raw_analyticdatasets_sas_participantstatusadult_pc"
#> [23] "cln_raw_analyticdatasets_sas_weeklevel_ee_pc"          
#> [24] "cln_raw_analyticdatasets_sas_weeklevel_re"             
#> [25] "cln_raw_biospec_collect_participant_assess"            
#> [26] "cln_raw_biospec_collect_participant_assess_24_hr"      
#> [27] "cln_raw_bld_pressure_heart_rate"                       
#> [28] "cln_raw_bld_spec_collect_tp1_pre_ex"                   
#> [29] "cln_raw_bld_spec_collect_tp2_20_min_ex"                
#> [30] "cln_raw_bld_spec_collect_tp3_40_min_ex"                
#> [31] "cln_raw_bld_spec_collect_tp4_10_min_post_ex"           
#> [32] "cln_raw_bld_spec_collect_tp5_30_min_post_ex"           
#> [33] "cln_raw_bld_spec_collect_tp6_3_5_hr_post_ex"           
#> [34] "cln_raw_bld_spec_collect_tp7_24_hr_post_ex"            
#> [35] "cln_raw_cardiopulmonary_ex_test"                       
#> [36] "cln_raw_ces_d"                                         
#> [37] "cln_raw_control_adherence_call"                        
#> [38] "cln_raw_control_monitoring_visit"                      
#> [39] "cln_raw_control_rest_record"                           
#> [40] "cln_raw_demographics"                                  
#> [41] "cln_raw_dhq_iii"                                       
#> [42] "cln_raw_dhq_iii_r"                                     
#> [43] "cln_raw_dxa"                                           
#> [44] "cln_raw_dxa_analysis"                                  
#> [45] "cln_raw_dxa_scan_results_ge"                           
#> [46] "cln_raw_dxa_scan_results_hologic"                      
#> [47] "cln_raw_dxa_scan_worksheet"                            
#> [48] "cln_raw_endurance_ex_tracking_log"                     
#> [49] "cln_raw_endurance_familiarization_session1"            
#> [50] "cln_raw_endurance_familiarization_session2"            
#> [51] "cln_raw_event_ascertainment"                           
#> [52] "cln_raw_grip_strength"                                 
#> [53] "cln_raw_height_weight_waist_circumference"             
#> [54] "cln_raw_intervention_monitoring_visit"                 
#> [55] "cln_raw_isometric_knee_extension"                      
#> [56] "cln_raw_local_lab_collect"                             
#> [57] "cln_raw_local_lab_results"                             
#> [58] "cln_raw_medical_history"                               
#> [59] "cln_raw_medication_inventory"                          
#> [60] "cln_raw_missed_ex_session"                             
#> [61] "cln_raw_mus_spec_collect_tp1_pre_ex"                   
#> [62] "cln_raw_mus_spec_collect_tp2_15_min_post_ex"           
#> [63] "cln_raw_mus_spec_collect_tp3_3_5_hr_post_ex"           
#> [64] "cln_raw_mus_spec_collect_tp4_24_hr_post_ex"            
#> [65] "cln_raw_participant_consent_status_log"                
#> [66] "cln_raw_participating_relatives"                       
#> [67] "cln_raw_pre_activity_medical_clearance"                
#> [68] "cln_raw_pre_screening_assess"                          
#> [69] "cln_raw_promis"                                        
#> [70] "cln_raw_randomization_enrollment_clearance"            
#> [71] "cln_raw_resistance_ex_tracking_log"                    
#> [72] "cln_raw_resistance_familiarization_session_1"          
#> [73] "cln_raw_resistance_familiarization_session_2"          
#> [74] "cln_raw_resistance_familiarization_session_3"          
#> [75] "cln_raw_resting_ecg"                                   
#> [76] "cln_raw_siteequipmentlog"
```

Even though you can directly call a data item, if you ever forget what
data items are available or want to know general categories of items,
use the load_clinical_data function.

``` r
clin_data = load_clinical_data()
names(clin_data)
#> [1] "curated"   "raw"       "chemistry"
names(clin_data$chemistry)
#> [1] "cln_chemistry_t02_plasma_lab_ck_results"  
#> [2] "cln_chemistry_t02_plasma_lab_conv_results"
#> [3] "cln_chemistry_t02_plasma_lab_crt_results" 
#> [4] "cln_chemistry_t02_plasma_lab_glc_results" 
#> [5] "cln_chemistry_t02_plasma_lab_ins_results"
```

The clinical files are separated into three different categories:

the curated files - representing the “main” datasets with the most
important and commonly used datasets

the raw files - representing the bulk of the data, most of these files
hold extremely detailed logistics of everything that went through
processing, storing, shipping, etc.

the chemistry files - representing clinical measurements for more
commonly measured chemical analytes. Includes glucose, insulin, creatine
kinase, etc etc.

### Omic modeling summary statistics (Differential Analysis)

See more documentation via `?load_differential_analysis`

``` r
differential_analysis = load_differential_analysis(
  repo_local_dir = NULL,
  selected_omes = "all",
  selected_tissues = "all",
  single_matrix = FALSE,
  epigen = FALSE,
  gsutil = "gsutil",
  combine_with_featgene = FALSE,
  verbose = TRUE
)
#> Please remember that the lowest CV Metabolite is chosen and the
#>             relevant refmet name is used. If you're not able to find your desired
#>             metabolite, look through the METABOLOMICS_CV object for the relevant
#>             refmet/feature name.
names(differential_analysis)
#> [1] "adipose" "blood"   "muscle"
names(differential_analysis[["blood"]])
#> [1] "metab"              "prot-ol"            "transcript-rna-seq"
```

By default,`load_differential_analysis` loads in the dataset in a nested
list in exactly the same way as the QC/expression datasets. Choose
whichever tissues or omes you’d like via `selected_omes` or
`selected_tissues`. You can find available tissues via
`tissue_available_list()` or `ome_available_list()`. Or if you enter in
a wrong mistaken tissue/ome, a warning or error will help.

If you would instead like to stack the matrixes more easily, use the
`single_matrix` function, which basically unlists the list and sticks
everything into a data.frame object.

``` r
single_matrix = load_differential_analysis(single_matrix = TRUE)
#> Please remember that the lowest CV Metabolite is chosen and the
#>             relevant refmet name is used. If you're not able to find your desired
#>             metabolite, look through the METABOLOMICS_CV object for the relevant
#>             refmet/feature name.
colnames(single_matrix)
#>  [1] "tissue"             "assay"              "platform"          
#>  [4] "full_model"         "contrast"           "contrast_short"    
#>  [7] "contrast_type"      "contrast_category"  "randomGroupCode"   
#> [10] "Timepoint"          "feature_id"         "logFC"             
#> [13] "CI.L"               "CI.R"               "degrees_of_freedom"
#> [16] "logLik"             "t"                  "AveExpr"           
#> [19] "z.std"              "p_value"            "adj_p_value"
```

For a quick explanation of each of the columns, you can find this via
`?load_differential_analysis`

Importantly, this loads the differential analysis for each of
comparisons mentioned in the methods, including the comparison between
the endurance or resistance group relative to time, fasting, biopsy,
etc. matched controls, comparison between the endurance and resistance
groups directly, and finally comparisons within group without a matched
control.

The majority of the analysis is done via exercise groups relative to the
controls (e.g. EE-CON, RE-CON) Make sure you filter to whichever
category you prefer before continuing with analysis.

``` r
single_matrix %>% dplyr::pull(contrast_category) %>% unique()
#> [1] EE-CON  RE-CON  EE-EE   RE-RE   EE-RE   CON-CON
#> Levels: EE-CON RE-CON EE-EE RE-RE EE-RE CON-CON
```

### Other items

Here’s a list of other items that come up relatively often. The
documentation of these is still a work in progress.

### Enrichment Results

``` r
colnames(CAMERA_RESULTS)
#>  [1] "tissue"         "assay"          "contrast_type"  "contrast"      
#>  [5] "contrast_short" "collection"     "database"       "set_id"        
#>  [9] "set"            "set_short"      "set_size"       "set_size_DB"   
#> [13] "size_ratio"     "direction"      "t"              "df"            
#> [17] "z.std"          "p_value"        "adj_p_value"
```

Quick summary: CAMERA-PR is a method of enrichment that incorporates all
features to generate a comparison of the test statistics between
in-pathway vs out-of-pathway test statistics to see if the statistics
within pathway are significant.

This file structure is more or less just an enrichment level match for
the comparisons described in the single-matrix differential analysis
results, where all tissues and assays are included in all the analysis.

### Feature to gene file:

``` r
head(HUMAN_FEATURE_TO_GENE)
#> Key: <assay, feature_id>
#>              assay               feature_id entrez_gene gene_symbol
#>             <fctr>                   <fctr>      <fctr>      <fctr>
#> 1: epigen-atac-seq chr1:100006105-100007013       23443     SLC35A3
#> 2: epigen-atac-seq chr1:100009408-100009608       23443     SLC35A3
#> 3: epigen-atac-seq   chr1:10001014-10001214      116362        RBP7
#> 4: epigen-atac-seq chr1:100010489-100010728       23443     SLC35A3
#> 5: epigen-atac-seq chr1:100021498-100021698       23443     SLC35A3
#> 6: epigen-atac-seq chr1:100024572-100024772       23443     SLC35A3
#>       ensembl_gene custom_annotation relationship_to_gene uniprot refmet_name
#>             <fctr>            <fctr>                <num>  <fctr>      <fctr>
#> 1: ENSG00000117620            Intron                    0    <NA>        <NA>
#> 2: ENSG00000117620              Exon                    0    <NA>        <NA>
#> 3: ENSG00000162444            Intron                    0    <NA>        <NA>
#> 4: ENSG00000117620            Intron                    0    <NA>        <NA>
#> 5: ENSG00000117620            Intron                    0    <NA>        <NA>
#> 6: ENSG00000117620            3' UTR                    0    <NA>        <NA>
#>    kegg_id flanking_sequence
#>     <fctr>            <fctr>
#> 1:    <NA>              <NA>
#> 2:    <NA>              <NA>
#> 3:    <NA>              <NA>
#> 4:    <NA>              <NA>
#> 5:    <NA>              <NA>
#> 6:    <NA>              <NA>
```

The feature-to-gene map links each feature tested in differential
analysis to a gene, using Ensembl version 105 (mapped to GENCODE 39) as
the gene identifier source. Proteomics feature IDs (UniProt IDs) were
mapped to gene symbols and Entrez IDs using UniProt’s mapping
files.Epigenomics features were mapped to the nearest gene using the
ChIPseeker::annotatePeak() function with Homo sapiens Ensembl release
105 gene annotations. Gene symbols, Entrez IDs, and Ensembl IDs were
assigned to features using biomaRt version 2.58.2 (Bioconductor 3.18).
This file links all of the features included in any ome/tissue in our
analysis. Use this to see how some levels of omic analysis (e.g. ATAC,
RNAseq) may link up in terms of ome names.
