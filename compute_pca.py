"""____________________________________________________________________________
Script Name:        04_pca.py
Description:        Compute PCA with all band's combinations.
Requirements:       GDAL version >= "3.1.4", Dask, scipy, numpy
Outputs:            High resolution image.
____________________________________________________________________________"""
# Parameters
# -----------------------------------------------------------------------------
mul_img = "./data/wv3_img.tif"
out_dir = "./PCA"
bands = ["C", "B", "G", "Y", "R", "RE", "N1"] # Ordered from first band to last band

# Main script
# -----------------------------------------------------------------------------
# Import packages
import os
# Add pca functions
import pca

# Create directories
def create_dir(dir_path: str, readme: str = None):
    """
    Create directories from path
    :dirPath: String with the new directory
    :readme: String to pass in README file

    **Important**: Dir is not created if it exists.
    """
    system_path = os.path.normpath(dir_path)
    if (not os.path.exists(system_path)):
        # Create directory
        os.makedirs(system_path)

    # Write README file
    if readme != None:
        file = open(os.path.join(system_path,'README.txt'), 'w')
        file.write(readme)
        file.close()

# Compute PCA
# -----------------------------------------------------------------------------
create_dir(out_dir,
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
