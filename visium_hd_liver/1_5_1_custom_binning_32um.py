import os
import harpy as sp
import spatialdata as sd
import spatialdata_io
import json

first_run = False

if first_run:
    sdata = spatialdata_io.visium_hd('/data/gent/vo/000/gvo00070/vsc43831/spatial_datasets/Visium_HD_Liver',
                                    #fullres_image_file="microscope_image/SCA002_HE_hires_transposed.tif",
                                    annotate_table_by_labels=True)
        
    # Write to zarr
    sdata.write('/data/gent/vo/000/gvo00070/vsc43831/spatial_datasets/Visium_HD_Liver/Visium_HD_Liver_zarr')
    del sdata
    
sdata = sd.read_zarr('/data/gent/vo/000/gvo00070/vsc43831/spatial_datasets/Visium_HD_Liver/Visium_HD_Liver_zarr')
print(sdata)

json_file = "/data/gent/vo/000/gvo00070/vsc43831/spatial_datasets/Visium_HD_Liver/binned_outputs/square_016um/spatial/scalefactors_json.json"

with open(json_file, 'r') as f:
    scalefactors = json.load(f)

# Print the loaded JSON data
print(scalefactors)
spot_diameter_32um = scalefactors['spot_diameter_fullres'] * 2

# take shape of grid equal to full h&e image
se = sdata["Visium_HD_Liver_full_image"]["scale0"]["image"]
sdata = sp.im.add_grid_labels_layer(sdata, shape = se.shape[1:], size=spot_diameter_32um, output_shapes_layer="Visium_HD_Liver_square_032um",
                                    output_labels_layer="Visium_HD_Liver_square_032um_labels", grid_type="square")
print("Added grid")
print(sdata)

sdata = sp.tb.bin_counts( 
    sdata,
    table_layer="square_002um",
    labels_layer="Visium_HD_Liver_square_032um_labels",
    output_layer="square_032um",
    append = False,
    overwrite = True
)

print("Binned counts")
print(sdata)