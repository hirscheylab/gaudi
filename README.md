![](gaudi_img.png)

> The package is named `gaudi`, inspired by the renowned Spanish architect Antoni Gaudí, who was famous for his intricate and colorful mosaics. Just as Gaudí pieced together countless fragments of tiles to create stunning, cohesive artworks, this package brings together various pieces of multi-omics data to form a comprehensive picture of complex biological systems. The concept of a mosaic beautifully parallels the integration of diverse data types, illustrating how individual fragments (genes, proteins, metabolites, etc.) come together to reveal the larger, more intricate patterns within biological research.

# GAUDI

The `gaudi` package provides a streamlined solution for the integration and analysis of complex multi-omics data. Leveraging the power of UMAP (Uniform Manifold Approximation and Projection), this package offers researchers an efficient and intuitive approach to uncover hidden patterns and relationships within diverse biological datasets. Designed for ease of use and compatibility with existing R-based bioinformatics workflows, `gaudi` is ideal for both novice and experienced users looking to dive deeper into the world of multi-omics research.

## System Requirements

- **Operating Systems**: Tested on macOS 12+, Ubuntu 20.04+, Windows 10+
- **Dependencies**: 
  - R version 4.2+
  - `umap` R package
  - `dplyr`, `ggplot2`, and other CRAN packages (specified in the package dependencies)
- **Hardware**: Standard desktop computer; higher RAM recommended for large datasets.

## Installation

To install the latest GitHub version:

```r
# install.packages("devtools")
devtools::install_github("hirscheylab/gaudi")
```

Typical install time: Less than 5 minutes on a normal desktop computer.

## Demo

Here's a comprehensive example demonstrating GAUDI's capabilities for multi-omics integration and analysis:

```r
library(gaudi)
library(dplyr)
library(survival)

# 1. Data Integration
result <- gaudi_integrate(
    data_list = list(
        expression = expression_matrix,
        methylation = methylation_matrix,
        mirna = mirna_matrix
    ),
    n_clusters = 5,
    min_pts = 5
)

# 2. Clinical Association Analysis
clinical_results <- gaudi_clinical_analysis(
    integrated_data = result,
    clinical_data = clinical_info,
    variables = c("age", "gender", "treatment_history")
)

# 3. Survival Analysis
survival_results <- gaudi_survival_analysis(
    integrated_data = result,
    survival_data = survival_info,
    time_column = "Survival",
    event_column = "Death"
)

# 4. Pathway Enrichment
enrichment_results <- gaudi_enrichment(
    integrated_data = result,
    database = "REACTOME",
    pval_threshold = 0.05
)

# 5. Visualization
plot_gaudi_clusters(result)
plot_gaudi_survival(survival_results)
plot_gaudi_enrichment(enrichment_results)
```

## Instructions for Use

1. **Data Preparation**
   ```r
   # Prepare your multi-omics data matrices
   # Each matrix should have samples as rows and features as columns
   data_list <- list(
       expression = expression_matrix,  # Gene expression data
       methylation = methylation_matrix,  # DNA methylation data
       protein = protein_matrix  # Protein abundance data
   )
   ```

2. **Basic Integration**
   ```r
   # Integrate multi-omics data with default parameters
   integrated_results <- gaudi_integrate(
       data = data_list,
       n_clusters = 5,  # Number of expected clusters
       min_pts = 5     # Minimum points for cluster formation
   )
   ```

3. **Advanced Analysis**
   ```r
   # Clinical associations
   clinical_associations <- gaudi_clinical_analysis(
       integrated_data = integrated_results,
       clinical_data = clinical_info,
       variables = c("age", "stage", "grade")
   )
   
   # Survival analysis
   survival_results <- gaudi_survival_analysis(
       integrated_data = integrated_results,
       survival_data = survival_info,
       time_column = "days_to_event",
       event_column = "event_status"
   )
   
   # Biological enrichment
   enrichment_results <- gaudi_enrichment(
       integrated_data = integrated_results,
       database_path = "path/to/annotation.gmt",
       pval_threshold = 0.05
   )
   ```

4. **Visualization**
   ```r
   # Basic clustering visualization
   plot_gaudi_clusters(
       integrated_results,
       color_by = "cluster",
       labels = TRUE
   )
   
   # Survival curves
   plot_gaudi_survival(
       survival_results,
       risk_table = TRUE,
       conf_int = TRUE
   )
   
   # Enrichment plots
   plot_gaudi_enrichment(
       enrichment_results,
       top_n = 20,
       show_pvalues = TRUE
   )
   ```

For detailed examples and benchmarking, refer to our [benchmarking repository](https://github.com/hirscheylab/umap_multiomics_integration).

## License

GAUDI is licensed under the GNU General Public License v3.0. See the [LICENSE file](https://github.com/hirscheylab/gaudi/blob/main/LICENSE.md) in the repository for details.

## Citation

```bibtex
@article {Castellano-Escuder2024.10.07.617035,
	author = {Castellano-Escuder, Pol and Zachman, Derek K and Han, Kevin and Hirschey, Matthew D},
	title = {Interpretable multi-omics integration with UMAP embeddings and density-based clustering},
	elocation-id = {2024.10.07.617035},
	year = {2024},
	doi = {10.1101/2024.10.07.617035},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2024/10/11/2024.10.07.617035},
	eprint = {https://www.biorxiv.org/content/early/2024/10/11/2024.10.07.617035.full.pdf},
	journal = {bioRxiv}
}
```
