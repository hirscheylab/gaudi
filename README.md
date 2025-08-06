
![](gaudi_img.png)

> The package is named `gaudi`, inspired by the renowned Spanish architect Antoni Gaudí, who was famous for his intricate and colorful mosaics. Just as Gaudí pieced together countless fragments of tiles to create stunning, cohesive artworks, this package brings together various pieces of multi-omics data to form a comprehensive picture of complex biological systems. The concept of a mosaic beautifully parallels the integration of diverse data types, illustrating how individual fragments (genes, proteins, metabolites, etc.) come together to reveal the larger, more intricate patterns within biological research.

# GAUDI

The `gaudi` package provides a streamlined solution for the integration and analysis of complex multi-omics data. Leveraging the power of UMAP (Uniform Manifold Approximation and Projection), this package offers researchers an efficient and intuitive approach to uncover hidden patterns and relationships within diverse biological datasets. Designed for ease of use and compatibility with existing R-based bioinformatics workflows, `gaudi` is ideal for both novice and experienced users looking to delve deeper into the world of multi-omics research.  

## Installation

To install the latest GitHub version:

``` r
# install.packages("devtools")
devtools::install_github("hirscheylab/gaudi")
```

## Parameter Settings and Dataset Adaptation

GAUDI provides comprehensive parameter customization while offering empirically-derived defaults that enable robust analysis across diverse multi-omics contexts:

### Default Parameters

- **UMAP Parameters**
  - `n_neighbors` = 15 (UMAP default)
  - `n_components` = 4 for initial embeddings, 2 for final visualization
  - `metric` = "euclidean"
  - `min_dist` = 0.01

- **HDBSCAN Parameters**
  - Minimum cluster size automatically computed as 3% of sample count (minimum of 2 samples)
  - Adapts to varying dataset sizes while preventing singleton clusters

- **Feature Selection Parameters (metagenes)**
  - XGBoost (default): lambda=0, eta=0.5, gamma=50, max_depth=10, subsample=0.95
  - Random Forest option available for broader feature inclusion

### Data Preprocessing

#### Data Preparation

Prepare your multi-omics data matrices. Each matrix should have samples as rows and features as columns:

```r
omics_list <- list(
    expression = expression_matrix,  # Gene expression data
    methylation = methylation_matrix,  # DNA methylation data
    protein = protein_matrix  # Protein abundance data
)
```
    
#### Missing Data

- Features with zero standard deviation are automatically removed
- We recommend imputing missing values before using GAUDI
- For large-scale missingness (>20%), consider removing affected samples/features

#### Batch Effects

- Use `combine_omics = TRUE` to enable automatic batch correction
- Employs `ComBat` from `sva` package to correct omic-type effects

```r
result <- gaudi(omics_list, combine_omics = TRUE)
```

### Adapting to New Datasets

#### Dataset Size Considerations

- For small datasets (<100 samples): Consider reducing UMAP `n_neighbors`
- For large datasets (>1000 samples): May benefit from increased `n_neighbors`
- HDBSCAN minimum cluster size automatically scales with dataset size

#### Feature Selection

- Use XGBoost (default) for selective feature identification
- Switch to Random Forest for broader feature inclusion:

```r
result <- gaudi(omics_list, method = "rf")
```

#### Validation

- Use silhouette scores to assess clustering quality
- Compare survival differences between clusters
- Perform enrichment analysis on identified features

```r
# Check clustering quality
print(result@silhouette_score)
```

_Like most machine learning approaches, while these parameters provide robust performance across typical multi-omics datasets, users are encouraged to optimize them for their specific biological context through standard cross-validation and benchmarking procedures, particularly when analyzing datasets with unique characteristics or investigating novel biological phenomena._

## Method Performance and Validation

For an in-depth understanding of the GAUDI (Group Aggregation via UMAP Data Integration) method's performance and its comparative analysis with other leading multi-omics integration techniques, we encourage users to explore our dedicated benchmarking repository. This repository contains detailed benchmarks across various datasets, including simulated data, TCGA cancer datasets, single-cell datasets, and DepMap multi-omics data, providing valuable insights into the effectiveness and versatility of the GAUDI approach. Access the comprehensive benchmarks and results [here](https://github.com/hirscheylab/umap_multiomics_integration).

## System Requirements

### Hardware requirements

- `GAUDI` package requires only a standard computer with enough RAM to support the in-memory operations.
- No specialized hardware required (GPU acceleration optional but not necessary)

### OS Requirements

This package has been tested on the following systems:

- macOS: Sonoma (14.6.1)

### Typical Install Time on a "Normal" Desktop Computer

- Less than 5 minutes for installation of the main package
- Up to 15 minutes if all dependencies need to be installed from scratch

### Expected Run Time

- Small-medium datasets: less than 30 seconds
- Benchmarking examples with TCGA data: 5-10 minutes
- Runtime scales with dataset size and dimensionality
  
## License

GAUDI is licensed under the GNU General Public License v3.0. See the [LICENSE file](https://github.com/hirscheylab/gaudi/blob/main/LICENSE.md) in the repository for details.

## Citation

``` bibtex
@article{castellano2025gaudi,
  title={GAUDI: interpretable multi-omics integration with UMAP embeddings and density-based clustering},
  author={Castellano-Escuder, Pol and Zachman, Derek K. and Han, Kevin and Hirschey, Matthew D.},
  journal={Nature Communications},
  volume={16},
  pages={5771},
  year={2025},
  publisher={Nature Publishing Group},
  doi={10.1038/s41467-025-60822-1},
  url={https://doi.org/10.1038/s41467-025-60822-1}
}
```

