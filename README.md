# 🔥 Automating Wildfire Hazard with Google Earth Engine
## University of Arizona - GIST 909

[Heather Robbins](https://github.com/hrobbins-gis)  

### Overview
This project developed an **automated workflow for wildfire hazard assessment** using **Google Earth Engine (GEE)**. Model validation was performed using **R** and **QGIS** software.  The workflow integrates remote sensing data collection, raster normalization, and susceptibility modeling to streamline wildfire hazard modeling.

The workflow was tested using the 2024 Line Fire in the San Bernardino National Forest.

By automating data processing and map generation, this project seeks to **reduce the time and effort** required to produce accurate, repeatable wildfire hazard maps — supporting better decision-making for resource managers and emergency planners.

---
### 🌎 Objectives
- Automate data collection and preprocessing using **Google Earth Engine**.
- Normalize and weight multiple risk factors in **Google Earth Engine** to compute a composite **Wildfire Hazard Index (WHI)**.
- Demonstrate a scalable open-source workflow that can be adapted for different regions.

---
## Quick Start Guide

- Open `scripts/gee/wildfire_hazard_model.js` in Google Earth Engine Code Editor
- Draw/input your study area polygon
- Update START and END dates
- Run the script to generate hazard outputs
- Optional:
- - Uncomment added map layers for visualization
  - Export raster outputs to Google Drive

---
### 🧩 Workflow Summary

1. **Data Acquisition (GEE)**
   - Vegetation health (USGS - Landsat 8 and 9)
   - Land cover (ESA - WorldCover)
   - Slope (SRTM - DEM)
   - Vapor pressure deficit (Daymet)
   - Precipitation (Daymet)
   - Wind speed (gridMET)

2. **Wildfire Hazard Index Modeling (GEE)**
   - Normalize each raster layer (0–1 scale)
   - Apply weights to reflect each factor’s relative influence
   - Combine layers to calculate a composite Wildfire Hazard Index (WHI)
   - Export results

---
### ✅ Model Validation
Model validation was conducted using R to assess agreement between wildfire hazard outputs and observed fire perimeter and indices (e.g., dNBR).

- Statistical tests: Wilcoxon test
- Accuracy assessment: Confusion matrix  
- Effect size metrics included  

Full results and interpretation are provided in the accompanying research paper. This repository includes scripts and sample outputs to support reproducibility.

---

### 🗂️ Repository Structure
```plaintext
wildfire-hazard-gee/
│
├── README.md              
├── LICENSE
├── .gitignore
│
├── docs/                  
│   ├── abstract.md
│   ├── research.pdf           
│   └── workflow-diagram.png
│
├── scripts/
│   ├── README.md               
│   ├── gee/               
│   │   ├── aoi.js
│   │   ├── wildfire_hazard_model.js
│   │   └── dnbr.js
│   │
│   ├── r/                 
│   │   └── model_validation.R
│
├── results/
│   ├── maps/
│   │   └── wildfire_hazard_map.jpg
│   ├── figures/
│   │   ├── boxplot.png
│   │   └── dnbr_scatter.png
│   └── outputs/
│       └── hazard.tif
```

### ⚙️ Dependencies
- **Google Earth Engine (JavaScript API)**
- **QGIS**
- - Version ≥ 3.30
- **R**
   - caret
   - dplyr
   - ggplot2
   - terra

### 🧾 License

This project is licensed under the MIT License.
You are free to use, modify, and distribute this work with attribution.


### 📬 Contact

For inquiries or collaboration, please contact [@hrobbins-gis](mailto:hrobbins@barstow.edu)
