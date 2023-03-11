This code repository is associated with **Similarities and Differences in Alzheimer’s Disease Comorbidities in Racialized Populations identified from Electronic Medical Records**

Overview: Code for Running Analyses, with guidance on how to define cohort and obtain diagnoses. Follow below after generating venv using requirements.txt

1. Defining the cohort 

    A. The following criteria must be met when identifying patients with Alzheimer's dementia in an EHR database:
        1. Patients must be aged 65+
        2. Patients must be diagnosed with at least one of the following ICD-10-CM codes:
            1. G30.1
            2. G30.8
            3. G30.9
    
    B. Use MatchIt for propensity score matching to identify matched controls
        1. 1:2 ratio used
        2. Controls were matched on estimated age, identified sex, death status, and identifiied race and ethnicity (R&E)
        3. First round of PS matching: R_notebooks/PS_Matching_1.R
            This results in: Demographics/AD_control_demographics.csv
        4. Demographics/AD_control_demographics.csv file is processed by Misc_Notebooks.ipynb in order to make
        separate demographic files for patients with AD and respective controls
            This results in: Demographics/ad_demographics.csv and Demographics/con_demographics.csv
        5. Second round of PS matching for patients with AD gets processed by R_Notebooks/PS_MI_RE_ad.Rmd (using 
        Demographics/ad_demographics.csv file)
            This results in: Demographics/RE_MI_ad_demo.csv
        6. Second round of PS matching for control patients gets processed by R_Notebooks/PS_MI_RE_con_v2.Rmd (using Demographics/con_demographics.csv file)
            This results in: Demographics/RE_MI_con_demo.csv
            
2. Obtaining diagnoses

    Patients' ICD-9-CM and ICD-10-CM diagnoses retrieved from the EHR database were then mapped to phecodes using SQL. Diagnoses pulled from the condition_occurrence table were limited to those with domain_id = 'Condition' and concept_class_id = 'Clinical Finding' for consistency with more standardized OMOP databases. Diagnoses were then mapped to phecodes using table prepared in phecodes/phecodes.ipynb, which consolidates two phecode tables (mapping ICD-9-CM to phecodes and mapping ICD-10-CM to phecodes). This results in the following csv files:
        A. ad_diagnoses_phecodes_omop.csv
        B. con_diagnoses_phecodes_omop.csv
        
3. Running Analyses - Run these notebooks in order. Note that some of these notebooks assume that the validation cohort files are already available.

    A. UMAP: Refer to 1_LowDimRepresent_RE.ipynb for analysis, which includes code for visualizing Fig 2A-F.
    B. Dunn's test for UMAP: R_Notebooks/Dunn_UMAP_RE.R
    C. Common phenotypes all patients with AD have across racialized populations: Misc_Notebooks/common_phenotypes.ipynb
    D. Phenotype Differential Analysis: Refer to 2_CompareADControls_RE.ipynb for analysis, which includes code for obtaining Supplementary Data File 1 as well as visualizing Fig 3.
    E. Network Analyses: Refer to 3_PhenotypeNetworks_RE.ipynb for analyses, which includes code for visualizing Supp Fig 5-6
    F. Dunn's test for Network Metrics: R_Notebooks/Dunn_Networks_RE.R
    
    
4. Figures, Tables, and Supplementary Data Files not created in analyses notebooks. Note that some of these notebooks assume that the validation cohort files are already available.

    A. Table 1: R_notebooks/Table_1.R
    B. Fig 4: Extra_Figure_SuppData_Notebooks/Fig_4_Manhattan_Plots.ipynb
    C. Fig 5: Extra_Figure_SuppData_Notebooks/Fig_5_Comparing_diagnoses_UCSF_UCDDP.ipynb
    D. Fig 6: Extra_Figure_SuppData_Notebooks/Fig_6_EF_Network_Metrics_Comp_UCSF_UCDDP.ipynb
    E. Supplementary Data Files 6 and 7: Extra_Figure_SuppData_Notebooks/SuppData_6_7.ipynb