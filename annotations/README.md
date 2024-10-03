## Description

There are 2 codes for annotating VISIUM spots with pathological annotations taken from **OMERO** and one code to annotate Xenium defined cells (or by another segementation method). 

**categorical_ann2SR.py** - is a code which uses OMERO annotations to allocate each VISIUM spot to one or several annotation ROIs. If VISIUM spot is splitted between 2 ROIs - the program will chose the one ROI which includes centroid of the spot. If the spot is inside "ROI_A" which is inside "ROI_B", the spot will be allocated to both ROIs without information about their relationships: "ROI_A;ROI_B". Output of this code is cswv file with barcodes of VISIUM spots and their corresponding annotation ROIs.

**ann2SR.py** - is more advanced version of **categorical_ann2SR.py** as it contains table with fractional distribution of VISIUM spot (which has a certain size) across all ROIs. Also there is a hierarchical information about the ROIs, their level, "children" (ROIs within a given ROI) and "parents" (all ROIs that contain given ROI). Spatial table of fractional spational distribution is saved in ```adata.obs``` and hierarchical information about ROIs is saved in ``` adata.uns[‘rois_hierarchy’] ```

Both of the functions use the same input (csv file) and output (output folder path)

**ann2Xenium.py** - is allocating cells to manually annotated regions in categorical manner (similar to **categorical_ann2SR.py**). There is no hierarchical information about ROIs, and if one ROI is inside another, then annotations will be added as  "ROI_A; ROI_B". As input it takes configuration file (for example see conf_Xenium.yaml), and as output produces csv files which are copy of segmentation files with additional column - 'annotation'.

**ann2Xenium_batch.py** - allows to run **ann2Xenium_v2.py** (which is functionally the same as **ann2Xenium.py**) for a batch of datasets. You will have to run it for master configuration file (see as an example **conf_Xenium_batch.yaml**). User can specify its OMERO credentials in this configuration file, however there is a possibility to enter your credentials after programs starts to avoid saving your password in yaml file. For this specify *None* for *omero_username* or *omero_password* flags.

**ann2Xenium_v2.py** - is basically the same as **ann2Xenium.py**, the only difference is how do you pass your OMERO credentials. in this version you will need to store your credentials in configuration file (so job can be runned non-interactively). Then make sure the configuration file is well protected (or delete you OMERO credentials from configuration file just after the job is done) 
TODO: change the code the way that it can read environmental parameters from bash, in this case credentials will not be stored anywhere in the file

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
or 
```
python ann2Xenium.py path/to/conf/file.yaml 
```
After running you will be asked for OMERO credentials

In case you have many regions with slightly different names or typos in their name, that belongs to one class (for example: "Basilar artery", "basilar artery", "Basilar Artery", "Basilar arteiry" etc), you can manually create a dictionary csv with two columns: "original" and "unified" that can be used to correct names, so all ROIs with some differences in the names will appear as one class. As example see **Dictionary_annotations_example.csv**. In this case you'll have to run the code as:
```
python ann2SR.py path/to/csv path/to/output/folder path/to/Dictionary_annotations.csv
```

## Reading output files
Output anndata files saved in *h5ad* format. To open it use:
```
import anndata
adata = anndata.read_h5ad(path_to_adata)
```
