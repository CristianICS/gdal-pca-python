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

*Note: The number of combinations for 7-bands image is high (99). Combinations
are cropped to 10, but this can be undone replacing* `combis_subset` *by*
`combis_dict` *in line 59 inside* `compute_pca.py` *file*.

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

Lastly the stats are computed in csv or json format.

```json
{
    "eigenvals": [2.96324133257161, 0.032160123770838975, 0.004598543657550729], 
    "eigenvals_%": [98.77471108572034, 1.0720041256946324, 0.1532847885850243], 
    "eigenvectors": [
        [-0.5752291616678924, -0.7718603913163038, 0.2708190316131709],
        [-0.5798495044419998, 0.15124346839158867, -0.8005622808172077],
        [-0.5769627056970156, 0.6175410509317776, 0.5345625189338048]
    ]
}
```

```csv
"stat_name","band_1","band_2","band_3"
"eigenvals",2.96324133257161,0.032160123770838975,0.004598543657550729
"eigenvals_%",98.77471108572034,1.0720041256946324,0.1532847885850243
"eigenvectors_pc1",-0.5752291616678924,-0.7718603913163038,0.2708190316131709
"eigenvectors_pc2",-0.5798495044419998,0.15124346839158867,-0.8005622808172077
"eigenvectors_pc3",-0.5769627056970156,0.6175410509317776,0.5345625189338048
```

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