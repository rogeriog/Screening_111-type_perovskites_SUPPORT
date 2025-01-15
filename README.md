# Supporting Files for "Screening 111-type Perovskites" Paper

This repository contains the data and scripts used to support the paper titled "Screening 111-type Perovskites". It includes:

*   **Processed data** from the Open Quantum Materials Database (OQMD)
*   **Machine learning (ML) predictions** for 111-type perovskite structures.
*   **Scripts** for data processing, filtering, and analysis.

## File Descriptions

### Data Files
*   **`OQMD_one_hot_encoded_data.tar.xz`**: A compressed archive of the one-hot encoded OQMD data that was used to train the coverage analysis model.

*   **`perovsk111_30kSamples_fullfeaturized_coverage.csv`**: This file contains the results of the machine learning predictions and the coverage scores. It includes the following columns:
    *   `structure`: The chemical formula of the perovskite structure.
    *   `gap_class`: The predicted value for the band gap class 'semiconductor', based on MODNet model trained on the OMEGA+ROSA featturizer.
    *   `band_gap`: The predicted band gap energy (in eV), based on MODNet model trained on the OMEGA+ROSA featturizer..
    *   `stability`: The predicted stability score, based on MODNet model trained on the OMEGA+ROSA featturizer.
    *   `coverage_score`: A score indicating how well the elemental combinations in this structure are represented in the OQMD database based on the `get_coverage_from_OQMD.py` script.
    *   `mlip_stability`: The predicted stability score, based on machine learned interatomic potentials, in our case CHGNet.
    
    *Example rows:*
    ```
    structure,gap_class,band_gap,stability,coverage_score,mlip_stability
    Cs2K4Ag3Bi(Cl2O)6,0.72,0.9941473,0.4423089851106121,0.2399999999999999,0.1043225371444671
    Cs2K4Ag3Mo(Cl2O)6,0.8,1.0516589,0.6043049747399143,0.2628571428571428,0.08490707282715244
    Cs2K4Ag3Pb(Cl2O)6,0.8800000000000001,1.1520798,0.3849695152404276,0.2285714285714285,0.09041165614926075
    ...
    ```
    
*   **`filtered_exploration_df.csv`**: This file contains a subset of the `perovsk111_30kSamples_fullfeaturized_coverage.csv` file, filtered to identify structures suitable for further exploration. It includes structures with:
    *   `stability <= 0.015`
    *   `mlip_stability <= 0.045`
    *   `band_gap <= 3.5`
    *   `gap_class >= 0.5`
    *   `coverage_score >= 0.25`

*   **`filtered_production_df.csv`**: This file contains a more strictly filtered subset of the `perovsk111_30kSamples_fullfeaturized_coverage.csv`  file, with structures considered for production. It includes structures with:
    *   `stability <= 0.015`
    *   `mlip_stability <= 0.045`
    *   `band_gap <= 3.5`
    *   `gap_class >= 0.5`
    *   `coverage_score >= 0.45`

### Python Scripts

*   **`get_coverage_from_OQMD.py`**: This script calculates the coverage score for each structure by analyzing the representation of elemental combinations within the OQMD.
    *   **Key Functionality:**
        *   Fetches data from a MySQL database containing OQMD entries.
        *   Removes polymorphs, keeping the most stable entry for each composition.
        *   Normalizes and reduces chemical compositions to unique sets of elements.
        *   Creates a one-hot encoded DataFrame of the OQMD data.
        *   Calculates a 'coverage score' for each composition in the target dataset by examining the coverage of element combinations.
        *   Outputs the results to the `perovsk111_30kSamples_fullfeaturized_coverage.csv` file.
        
    *   **Algorithm Highlights:**
        *   It utilizes a lookup table for the number of possible elemental combinations for 2 to 8 element structures.
        *   The coverage score is calculated by checking the presence of these combinations in the OQMD dataset, giving more weight to higher combinations of elements in the coverage calculation.
        *   The script supports incremental processing of structures to prevent memory overload, and stores the progress in `progress.json`.

*   **`filter_ML_screen_predictions.py`**: This script filters the machine learning predictions based on several criteria, generating the `filtered_exploration_df.csv` and `filtered_production_df.csv` files.
    *   **Key Functionality:**
        *   Reads the `perovsk111_30kSamples_fullfeaturized_coverage.csv` file.
        *   Applies filters based on `stability`, `mlip_stability`, `band_gap`, `gap_class`, and `coverage_score`.
        *   Creates two filtered dataframes: `filtered_exploration_df` and `filtered_production_df`, each for different purposes.
        *   Groups the `filtered_production_df` by the B-site elements and outputs this information.

## Usage

1.  **Data Preparation:** Ensure you have access to the OQMD database if you want to repoduce the element encoded dataset obtained with `get_coverage_from_OQMD.py`.
2.  **Coverage Calculation:** Run `get_coverage_from_OQMD.py` to calculate the coverage scores of the structures in `CSIp3m1_30kselected_relaxed_featurized_OmegaROSA_predictions.csv`
3.  **Filtering:** Run `filter_ML_screen_predictions.py` to filter the structures into exploration and production datasets.

## Additional Notes
*   The  `perovsk111_30kSamples_fullfeaturized_coverage.csv` file contains all structures that were screened by the active learning model using the protostructure ML model along with the
    decomposition energies calculated through CHGNet and after CHGNet relaxation. Structures with Echgnet_stab above 45 meV are included in this set but they are not part of the screening,
    they just helped the model generalization, the subsets that matter are in `filtered_exploration_df.csv` and `filtered_production_df.csv`.
*   The `OQMD_one_hot_encoded_data.tar.xz` file is needed to generate new coverage scores based on OQMD through `get_coverage_from_OQMD.py`.
*   The parameters in the `filter_ML_screen_predictions.py` script can be adjusted to modify the filtering criteria based on the required exploration/production needs.
