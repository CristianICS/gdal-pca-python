"""____________________________________________________________________________
Script Name:          warp_pca.py
Description:          Compute gdalwarp tool in Osgeo shell to PCA rasters. 
                      GDAL version "3.1.4"
Outputs:              Clipped PCA rasters.
____________________________________________________________________________"""
import os

# Directorio principal
input_path = os.path.join("E:\\arqueologia_2022","celsa_dron",
"Yacimiento_Lepida_Celsa","Lepida Celsa","Mapas de reflectancia")

# Funcion que devuelve los parametros necesarios en cada zona
def setZoneParameters(input_path, name):

    # Directorio con las bandas
    input_zone_path = os.path.join(input_path, name, 'ACP')

    # Directorio de salida
    output_path_clipped = os.path.join(input_path, name, 'ACP_clip')

    # Prefijo de los archivos raster con las bandas
    prefix_bandname = "Lepida_Celsa_" + name + "_ACP"

    # Capa del geojson a utilizar como mascara (dependiendo de la zona)
    mask_layer = name.lower()

    return(
        input_zone_path,
        output_path_clipped,
        prefix_bandname,
        mask_layer
    )

# IDs de las combinaciones (1 - 42)
ids = range(1,43)

zones = ["Yacimiento"]

# gdalwarp tool
# ------------------
# Igualar el área de todas las bandas
for zone_list in zones:
    print("\nClip data\n=======================\n")

    # Save zone name
    zone = zone_list
    print("Zona "+zone+"\n")
    # Set zone parameters
    ip_zone, op_clip, prefix, mask_lyr = setZoneParameters(input_path, zone)

    # Archivo GPK con las mascaras para realizar gdalwarp
    mask_path = os.path.join(input_path, "masks.gpkg")

    # GDAL command
    command = (
        'gdalwarp -of GTiff' +
        ' -cutline "' + mask_path +
        '" -cl ' + mask_lyr + ' -crop_to_cutline' +
        ' -co COMPRESS=LZW -co PREDICTOR=3'+
        ' "{input}" "{output}"'
    )

    for file in os.listdir(ip_zone):
        if file.endswith('.tif'):
            gdal_input = os.path.join(ip_zone, file)
            filename = os.path.splitext(os.path.basename(file))[0]
            print("\nClip file "+filename+"\n")
            gdal_output = os.path.join(op_clip, filename + '.tif')
            # Run GDAL command in terminal
            os.system(command.format(input=gdal_input, output=gdal_output))