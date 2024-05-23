# Chromatography Technique Selector

## Project Overview

This project automates the selection of the most suitable chromatography technique for analyzing a given mixture of organic compounds. By simply providing the names of the molecules in the mixture, the code retrieves their physicochemical properties from various web sources and determines the optimal chromatography method based on these properties.

## How It Works

1. **Input**: 
   - User provides the names of the molecules present in the mixture.

2. **Data Retrieval**:
   - The code searches online databases and web pages to gather key physicochemical properties for each molecule:
     - **Boiling temperature (°C)**
     - **logP (partition coefficient)**
     - **pKa (acid dissociation constants)**
     - **Molecular mass**

3. **Chromatography Type Decision**:
   - **Gas Chromatography (GC)**:
     - Chosen if the minimum boiling temperature of the compounds is 300°C or higher.
   - **Ion Chromatography (IC)**:
     - Selected if the maximum molecular mass is less than or equal to 2000, and the maximum logP is negative, with a proposed pH derived from the pKa values.
   - **High-Performance Liquid Chromatography (HPLC)**:
     - Applied under various conditions depending on logP values and pKa ranges, with different stationary phases and eluent natures suggested based on the specific properties.
   - **Size Exclusion Chromatography (SEC)**:
     - Used for compounds with a molecular mass greater than 2000. The eluent nature and exact method depend on the logP values:
       - For positive logP values, SEC on gel permeation with an organic solvent is selected.
       - Otherwise, SEC on gel filtration with a polyhydroxylated hydrophile polymer in an aqueous solvent is recommended.

4. **Output**:
   - The code outputs the advisable chromatography type, the nature of the eluent (gas, aqueous, or organic), and the proposed pH for the eluent if applicable.

## Benefits

- **Efficiency**: Automates the property retrieval and decision-making process, saving time and reducing manual effort.
- **Accuracy**: Utilizes precise physicochemical data to ensure the most suitable chromatography technique is chosen.
- **Versatility**: Supports a wide range of organic compounds and chromatography methods.

## Getting Started


## Contributing
Contributions are welcome! Please submit a pull request or open an issue to discuss any changes.



