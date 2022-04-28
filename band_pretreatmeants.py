"""____________________________________________________________________________
Script Name:          band_pretreatmeants.py
Description:          Compute gdalwarp tool in Osgeo shell. GDAL version "3.1.4"
Outputs:              Resized & Clipped rasters.
____________________________________________________________________________"""
import os

# Directorio principal
input_path = os.path.join("E:\\arqueologia_2022","celsa_dron",
"Yacimiento_Lepida_Celsa","Lepida Celsa","Mapas de reflectancia")

# Zones data
zones = [
    {
        'name': "Casa_Hercules",
        # VISNIR bands parameters
        'resolution': 0.05163000000000000228,
        'bounds': '714514.8251200000522658 4583187.7234300002455711 714728.3668000000761822 4583374.9854399999603629'
    },
    {
        'name': "Yacimiento",
        # Thermal band parameters
        'resolution': 0.119990000000000,
        'bounds': '714341.3273200000403449 4582740.8748300001025200 715477.1526600000215694 4583983.0113099999725819'
    },
    {
        'name': "Campos",
        # Thermal band parameters
        'resolution': 0.119950000000000,
        'bounds': '713994.8637900000903755 4580924.6229500006884336 714980.3729900000616908 4583208.2310500005260110'
    }
]

# Funcion que devuelve los parametros necesarios en cada zona
def setZoneParameters(input_path, name):

    # Directorio con las bandas
    input_zone_path = os.path.join(input_path, name)

    # Directorio de salida
    output_path_clipped = os.path.join(input_path, name, 'clipped')
    # Directorio de salida de las capas redimensionadas
    output_path_resized = os.path.join(input_path, name, 'resized')

    # Prefijo de los archivos raster con las bandas
    prefix_bandname = "Lepida_Celsa_" + name + "_"

    # Capa del geojson a utilizar como mascara (dependiendo de la zona)
    mask_layer = name.lower()

    return(
        input_zone_path,
        output_path_clipped,
        output_path_resized,
        prefix_bandname,
        mask_layer
    )

# Bandas
bandas = [
    ["ndvi", "MS_index_ndvi"],
    ["green", "MS_transparent_reflectance_green"],
    ["red", "MS_transparent_reflectance_red"],
    ["red_edge", "MS_transparent_reflectance_red_edge"],
    ["nir", "MS_transparent_reflectance_nir"],
    ["thermal", "transparent_reflectance_thermal_ir_withoutOutlayers"]
]

# Localizar los scripts de GDAL dentro de Osgeo4W Shell
# os.system("py3_env")

# gdalwarp tool
# ------------------
# Igualar resolucion y tamaño de todas las bandas
for zone_list in zones:
    print("\nResize data\n======================\n")
    # Save zone name
    zone = zone_list['name']
    print("Zona "+zone+"\n")
    # Set zone parameters
    ip_zone, op_clip, op_res, prefix, mask_lyr = setZoneParameters(input_path, zone)

    # GDAL command
    res = zone_list['resolution']
    command = (
        'gdalwarp -of GTiff -te ' + zone_list['bounds'] +
        ' -tr ' + str(res) + ' -' + str(res) +
        ' -r cubicspline' +
        ' -co COMPRESS=LZW -co PREDICTOR=3'+
        ' "{input}" "{output}"'

    )

    # Apply above command in tif elements (UAV bands)
    for banda in bandas:
        print("\nResize "+banda[0]+"\n")
        gdal_input = os.path.join(ip_zone, prefix + banda[1] + '.tif')
        gdal_output = os.path.join(op_res, prefix + banda[1] + '_resize.tif')

        # Run GDAL command
        os.system(command.format(input=gdal_input, output=gdal_output))

# gdalwarp tool
# ------------------
# Igualar el área de todas las bandas
for zone_list in zones:
    print("\nClip data\n=======================\n")

    # Save zone name
    zone = zone_list['name']
    print("Zona "+zone+"\n")
    # Set zone parameters
    ip_zone, op_clip, op_res, prefix, mask_lyr = setZoneParameters(input_path, zone)

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

    for file in os.listdir(op_res):
        if file.endswith('.tif'):
            gdal_input = os.path.join(op_res, file)
            filename = os.path.splitext(os.path.basename(file))[0]
            print("\nClip file "+filename+"\n")
            gdal_output = os.path.join(op_clip, filename + '_clip.tif')
            # Run GDAL command in terminal
            os.system(command.format(input=gdal_input, output=gdal_output))
