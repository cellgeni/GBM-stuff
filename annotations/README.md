## Description

There are 2 codes for annotating VISIUM spots with pathological annotations taken from **OMERO**. 

**categorical_ann2SR.py** - is a code which uses OMERO annotations to allocate each VISIUM spot to one or several annotation ROIs. If VISIUM spot is splitted between 2 ROIs - the program will chose the one ROI which includes centroid of the spot. If the spot is inside "ROI_A" which is inside "ROI_B", the spot will be allocated to both ROIs without information about their relationships: "ROI_A;ROI_B". Output of this code is cswv file with barcodes of VISIUM spots and their corresponding annotation ROIs.

**ann2SR.py** - is more advanced version of **categorical_ann2SR.py** as it contains table with fractional distribution of VISIUM spot (which has a certain size) across all ROIs. Also there is a hierarchical information about the ROIs, their level, "children" (ROIs within a given ROI) and "parents" (all ROIs that contain given ROI). Spatial table of fractional spational distribution is saved in ```adata.obs``` and hierarchical information about ROIs is saved in ``` adata.uns[‘rois_hierarchy’] ```

Both of the functions use the same input (csv file) and output (output folder path)



## Preparation

Install conda environment using **environment.yml**. Make sure you have admin credentials for OMERO before running the script. Make a csv table file with information about OMERO ID of the sample, path to SpaceRanger output and rotational information. See **example.csv**

To install, run:
```
conda env create -n your_env_name -f environment.yml
```

## Running

Activate conda environment. Run python script as:
```
python ann2SR.py path/to/csv path/to/output/folder 
```
or:
```
python categorical_ann2SR.py path/to/csv path/to/output/folder 
```

After running you will be asked for OMERO credentials
