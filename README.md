# Bumblebee wing FA pipeline (R)

This repo contains the workflow for the analysis of **forewing fluctuating asymmetry (FA)** based on **StereoMorph** and **geomorph**.  

## What this pipeline does
- Processes wing photos and standardizes orientation
- Digitizes landmarks in StereoMorph
- Adds QC checks for the digitized landmarks
- Builds geomorph coordinate arrays + metadata tables
- Extracts FA (using single digitization per wing)
- Estimates digitization error from replicate digitizations
- Fits a simple FA ~ urbanity mixed model

## Folder structure (expected)

```text
images/
  raw_images/
  processed_forewings/
    lapidarius/
    pascuorum/
  processed_hindwings/                 # optional/legacy
    lapidarius/
    pascuorum/

data/
  landmarks_forewings.txt              
  landmarks_hindwings.txt              # optional/legacy
  landmark_data_forewings/
    lapidarius/
    pascuorum/
  landmark_data_hindwings/             # optional/legacy
    lapidarius/
    pascuorum/
  fw_coords_<species>.RDS              
  fw_metadata_<species>.csv            
  hw_coords_<species>.RDS              # optional/legacy
  hw_metadata_<species>.csv            # optional/legacy

results/
  05_fa_analysis/
  06_landmark_sensitivity/
  07_digitization_error/
  08_fa_models/
```
---

## Run order (scripts)

1. `01_process_images.R` — process + sort images  
2. `02_digitize_wings.R` — digitize landmarks in StereoMorph
3. `03_landmarks_qc.R` — QC: landmark count + likely swap flags  
4. `04_build_arrays.R` — build geomorph arrays + metadata
5. `05_extract_fa.R` — extract FA 
6. `06_landmark_sensitivity.R` — compare FA using 21 vs 15 landmarks (optional)  
7. `07_digitization_error.R` — digitization error from series 1 vs 2 (optional)  
8. `08_fa_urbanity_model.R` — FA ~ urbanity mixed model  

