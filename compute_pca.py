"""____________________________________________________________________________
Script Name:        04_pca.py
Description:        Compute PCA with all band's combinations.
Requirements:       GDAL version >= "3.1.4", Dask, scipy, numpy
Outputs:            High resolution image.
____________________________________________________________________________"""
# Parameters
# -----------------------------------------------------------------------------
mul_img = "D:/arqueologia_2022/wv3/015105921010_01/17MAY29064346/Pansharpen/17MAY29064346-M2AS-015105921010_01_P001_wBrovey.tif"
out_dir = "D:/arqueologia_2022/wv3/015105921010_01/17MAY29064346/PCA"
bands = ["C", "B", "G", "Y", "R", "RE", "N1"] # Ordered from first band to last band

# Main script
# -----------------------------------------------------------------------------
# Import packages
import os, sys
# Add modules
modules_path = os.path.abspath(os.path.join(os.path.normpath(__file__), os.pardir, 'modules'))
sys.path.append(modules_path)
import handlefiles as hf
import pca

# Compute PCA
# -----------------------------------------------------------------------------
hf.create_dir(out_dir,
    "Images with Principal Components.\n" +
    "Image is named with the combination id which identifies the bands with PCA" +
    "has computed.\n" +
    "Each image has the three componentes with more explained variancia. The first" +
    "component is located in the first band and so on.\n"
)

combis_dict = pca.combis(bands, 3, out_dir)

for combi_id, combi_bands in combis_dict.items():
    # Create PCA image
    print(f"\nWriting PCA from combi {combi_id} bands\n")
    print(f"==========================================\n")
    img_out_name = os.path.basename(mul_img)[0:39] + '_PCA_' + str(combi_id) + '.tif'
    img_out_path = os.path.join(out_dir, img_out_name)
    eigenvals, eigenvectors = pca.compute_pca(mul_img, combi_bands, img_out_path)
    # Write stats
    out_stats = os.path.join(out_dir, 'stats_combi_'+ str(combi_id) +'.csv')
    pca.write_stats(eigenvals, eigenvectors, out_stats, "csv")
