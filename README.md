# PCA with GDAL & Python

Apply Principal Component Analysis (PCA) over image bands. The function to
compute PCA is applied on all posible bands combination without repetition.
Combinations are computed with `itertools` for a minimun of 3 bands and
maximun same as total number of bands.

The process is performed by `computed_pca.py`; open this file and write the
following parameters:

- `mul_img` = Path to multispectral image.
- `out_dir` = Output directory folder.
- `bands`   = Image custom band keys in the same order as image bands.

Note: `pca.py` contains the PCA required functions.

## Procedure

First, the output directory is created and the band combinations are computed.

PCA is derived for each band combination and a new image, shaped with
the three first components, is exported. *Note: band combinations are
referenced with an ID, which is used to named the new PCA image.*

The main function is based on
[PCA4CD QGIS plugin](https://github.com/SMByC/PCA4CD/blob/master/core/pca_dask_gdal.py),
created by Xavier Corredor Llano.

1. Transform required image bands into a numpy array. Each band is flatten to
convert the data into a dataframe structure.
2. Group each flattened bands into a new array.
3. Compute covariance/correlation matrix.
4. Get *eigenvectors* and *eigenvalues* with `np.linalg.eigh()` and the matrix
from step 3.
5. Rearrange *eigenvalues* and *eigenvectors* started from PC with more
variance.
6. Apply *eigenvectors* to the data. This rotates and scales the data.
The principal components are now aligned with the axes of image bands.
7. Export the first three principal components.

## Workflow (Recommended)

Install [miniconda](https://docs.conda.io/en/latest/miniconda.html) and
create new environment:

```
conda create -n pca
```

Then activating the environment with `conda activate pca`, and install
following packages:

- GDAL
    ```
    conda install -c conda-forge gdal
    ```
- Dask
    ```
    conda install -c conda-forge dask
    ```

Finally run code with `python path/to/compute_pca.py`.

## Test data

WorldView3 clipped image. It's derived from a pansharpened product
(0.3m pixel size) computed with "weighted" Brovey method. It contains 
seven bands: "C", "B", "G", "Y", "R", "RE", "N1".

The image has radiometric and atmospheric corrections (DOS method - Chavez, 1989).

![Example of original image and PCA derived product.](data-example.png)

## Further development

- Scaling the data.

## References

- [A Step-by-Step Explanation of Principal Component Analysis (PCA)](https://builtin.com/data-science/step-step-explanation-principal-component-analysis)
-[The Basics: Principal Component Analysis](https://towardsdatascience.com/the-basics-principal-component-analysis-83c270f1a73c?gi=84c269d8c697)
- [PCA with R](https://www.datacamp.com/community/tutorials/pca-analysis-r)
- [Does mean centering or feature scaling affect a Principal Component Analysis?](https://sebastianraschka.com/faq/docs/pca-scaling.html)
- [Eigenvalues and eigenvectors in PCA](https://towardsdatascience.com/eigenvalues-and-eigenvectors-378e851bf372)