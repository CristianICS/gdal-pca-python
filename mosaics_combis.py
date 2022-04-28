"""____________________________________________________________________________
Script Name:          mosaics_combis.py
Description:          Calculate all posible combinations to create mosaics with
                      the given bands. GDAL version "3.1.4"
Outputs:              CSV file with all combis. 
____________________________________________________________________________"""
import os
import itertools

input_path = os.path.join("E:\\arqueologia_2022","celsa_dron","scripts","GDAL_R")

bandas = {
    "ndvi": "MS_index_ndvi",
    "green": "MS_transparent_reflectance_green",
    "red": "MS_transparent_reflectance_red",
    "red_edge": "MS_transparent_reflectance_red_edge",
    "nir": "MS_transparent_reflectance_nir",
    "thermal": "transparent_reflectance_thermal_ir_withoutOutlayers"
}

# Check if there is a thermal band
def checkThermal(keys):
    treatmeant = ''
    if("thermal" in keys):
        # Use clipped bands
        treatmeant += 'Resize-Compress-Clip'
    else:
        # Use resized bands
        treatmeant += 'None'
    
    return(treatmeant)

# Write metadata: bands stored in the mosaic + mosaic's id
def writeCombi(id, keys):
    # Bands inside the combination
    combi_bands = ' '.join(keys)

    # Treatmeant info (if there is no thermal band,
    # catch compressed original bands)
    band_treatmeant = checkThermal(keys)

    # Write line to insert in mosacis_ids file
    line = '\n'+str(id)+',"'+combi_bands+'","'+band_treatmeant+'"'

    return(line)

# Numero de bandas a combinar
n_combi = list(range(3,8))
# ID que identifique cada combinacion
combi_id = 1

# Open CSV files
csvOutput = open(os.path.join(input_path,"mosaics_ids.csv"),"w")
# Write header
header = '"id","band_name","band_pretreatmeant"'
csvOutput.write(header)

for L in n_combi:
    for subset in itertools.combinations(list(bandas.keys()), L):
        
        # Subset = list with subset's keys from bandas dictionary
        # bandas_subset = {key: bandas[key] for key in subset}

        combi = writeCombi(combi_id, subset)
        # Write the metadata of new raster
        csvOutput.writelines(combi)
        # Update id
        combi_id += 1

# Close csv 
csvOutput.close()

"""
for zone in zones:



    # Directorio con las bandas
    input_zone_path = os.path.join(input_path, zone_suffix)
    # Directorio de salida
    output_path = os.path.join(input_path, zone, 'mosaics')

    # Prefijos de los archivos raster con las bandas
    prefix_bandname = "Lepida_Celsa_" + zone + "_"

    # Numero de bandas a combinar
    n_combi = list(range(3,8))

    # ID que identifique cada combinacion
    id = 1

    # Txt with mosaics bands combinations
    #legend = open(os.path.join(output_path, 'legend.txt'), 'w')

    legend = open('legend.txt', 'w')
    



# Lineas de codigo a integrar en el Shell por numero de bandas
    # header = ['\ngdalwarp - ' + zone, '\n', '_'*50]
    # txt_output.writelines(header)


# Return list with the selected band's names
def computeRasters(band_names, prefix, suffix, input_path):
    rasters = []

    for i in band_names.keys():
        new_raster = os.path.join(input_path, prefix + band_names[i][1] + suffix)
        rasters.append(new_raster)

    return (rasters)

# Write metadata: bands stored in the mosaic + mosaic's id
def writeLegend(id, subset, bands):
    # Write line to insert in legend.txt
    line = id + ": "

    for i in subset:
        line = line + bands[i][0] + "  "

    # Add enter
    line = line + "\n"
    return line

# Loop to:
# 1. Generate mosaic combinations for each zone
# 2. Create mosaics with all combis and lengths
#
for zone in zones:

    print("Procesar bandas de "+zone+"\n__________________________")

    # Las bandas de 'Campos' no se han recortado
    if (zone == 'Campos'):
        zone_suffix = zone
        band_suffix = ".tif"
    else:
        zone_suffix = os.path.join(zone,"masked_bands")
        band_suffix = "_clip.tif"

    # Directorio con las bandas
    input_zone_path = os.path.join(input_path, zone_suffix)
    # Directorio de salida
    output_path = os.path.join(input_path, zone, 'mosaics')

    # Prefijos de los archivos raster con las bandas
    prefix_bandname = "Lepida_Celsa_" + zone + "_"

    # Numero de bandas a combinar
    n_combi = list(range(3,8))

    # ID que identifique cada combinacion
    id = 1

    # Txt with mosaics bands combinations
    #legend = open(os.path.join(output_path, 'legend.txt'), 'w')

    legend = open('legend.txt', 'w')
    for L in n_combi:
        for subset in itertools.combinations(list(bandas.keys()), L):
            # Subset = list with combination subset with band's keys
            bandas_subset = {key: bandas[key] for key in subset}

            # Lista con los directorios de las bandas
            rasters = computeRasters(bandas_subset, prefix_bandname, band_suffix, input_zone_path)

            # Write the metadata of new raster
            legend.writelines(writeLegend(str(id), subset, bandas))

            id += 1

    # Abrir raster
    rasterObjs, arrs, mask, projection, pixelResolution, intersectedBbox = raster.read(
        rasters)

    # pca transform
    eigValues, eigVectors, PCA = pca.pcaFnc(arrs)

    # save pca bands
    outfile = os.path.join(output_path, "pca.tif")
    if os.path.exists(outfile):
        os.remove(outfile)
    raster.write(PCA, mask, projection, pixelResolution,
                intersectedBbox, np.finfo(arrs.dtype).min, outfile)


    # Generar combinaciones con cada numero de bandas
    for L in n_combi:
        for subset in itertools.combinations(list(bands.keys()), L):
            # Output band
            output_band = os.path.join(output_path, prefix_bandname + "mosaic_" + str(id) + ".tif")
            # Generate GDAL command
            command = 'gdal_merge -separate -ot Float32 -of GTiff -co COMPRESS=DEFLATE -o "' + output_band + '" '
            # Subset = list with combination subset with band's keys
            command = createGDALMergeCommand(subset, command, input_zone_path, prefix_bandname, band_suffix)

            print(command+"\n\n")

            # Run in terminal
            os.system(command)

            # Write the metadata of new raster
            legend.writelines(writeLegend(str(id), subset, bands))

            print("Mosaico ID "+str(id)+" generado...\n")

            id += 1

rasterObjs, arrs, mask, projection, pixelResolution, intersectedBbox = raster.read(
    rasters)

# pca transform
values, vectors, projected = pca.pcaFnc(arrs)

# min-max normalization for projected array
for bidx in range(projected.shape[0]):
    min = np.min(projected[bidx])
    max = np.max(projected[bidx])
    projected[bidx] = (projected[bidx] - min) / (max - min)

# save pca bands
outfile = '/home/kikat/pca.tif'
if os.path.exists(outfile):
    os.remove(outfile)
raster.write(projected, mask, projection, pixelResolution,
             intersectedBbox, np.finfo(arrs.dtype).min, outfile)
"""