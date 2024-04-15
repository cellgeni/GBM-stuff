import os
import io
import json
import logging
from getpass import getpass

import omero
import numpy as np
import scanpy as sc
import pandas as pd
from tqdm import tqdm
from PIL import Image
from omero.gateway import BlitzGateway
from shapely import affinity,plotting
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
import matplotlib.pyplot as plt
import squidpy as sq
import fire

log = logging.getLogger(__name__)

def get_OMERO_credentials():
    logging.basicConfig(format='%(asctime)s %(message)s')    
    log.setLevel(logging.DEBUG)
    omero_host = "wsi-omero-prod-01.internal.sanger.ac.uk"
    omero_username = input("Username")
    omero_password = getpass("Password")
    return omero_host, omero_username, omero_password

def collect_ROIs_from_OMERO(omero_username, omero_password, omero_host, omero_image_id):
    ROIs = []
    log.info(f"Connecting to OMERO at {omero_host}")
    with BlitzGateway(omero_username, omero_password, host=omero_host, secure=True) as conn:
        # search image
        log.info(f"Looking for ImageId {omero_image_id}")
        conn.SERVICE_OPTS.setOmeroGroup('-1')
        image = conn.getObject("Image", omero_image_id)
        # set group for image
        log.info(f"Found image id={image.id} name='{image.name}'")
        log.info(f"Found image in group id={image.details.group.id.val} name='{image.details.group.name.val}'")
        
        log.info(f"Storing rendered thumbnail in memory for QC")
        img_data = image.getThumbnail() #tiny preview image
        rendered_thumb = Image.open(io.BytesIO(img_data))
        
        group_id = image.details.group.id #check group id
        conn.setGroupForSession(group_id.val)
        # get image ROIs
        roi_service = conn.getRoiService()
        log.info("Retrieving ROIs")
        result = roi_service.findByImage(image.id, None)
        #result has property roi - for one given image id
        for roi in result.rois:
            try:
                primary_shape = roi.getPrimaryShape()
                name = primary_shape.getTextValue().val
                #, separates x and y and  space separates points
                points = [(lambda xy : list(map(float,xy.split(","))))(xy) for xy in primary_shape.getPoints().val.split(" ")]
                ROIs.append({
                    "name": name,
                    "points": points
                })
                log.debug(f"Found ROI id={roi.id.val} name='{name}' type={primary_shape.__class__.__name__}")
            except:
                pass
    log.info(f"Found {len(ROIs)} ROIs in total")
    return ROIs, image

def read_SR_to_anndata(spaceranger_path):
    log.info(f"Reading spaceranger output from: '{spaceranger_path}'.")
    adata = sq.read.visium(spaceranger_path)
    return adata

def rotate_flip_polygon(polygon, img_center, rot_angle, flipX=False, flipY=False):
    angle = np.deg2rad(rot_angle)
    move1_matrix = np.array([[1,0,img_center[0]],[0,1,img_center[1]], [0,0,1]])
    rot_matrix = np.array([[np.cos(angle), -np.sin(angle), 0],
                          [np.sin(angle), np.cos(angle), 0],
                          [0,0,1]])
    flip_matrix = np.eye(3,3)
    if flipX: flip_matrix[0,0] = -1
    if flipY: flip_matrix[1,1] = -1  
    if np.abs(rot_angle) == 90 or np.abs(rot_angle) == 270:
        move2_matrix = np.array([[1,0,-img_center[1]],[0,1,-img_center[0]], [0,0,1]])
    else:
        move2_matrix = np.array([[1,0,-img_center[0]],[0,1,-img_center[1]], [0,0,1]])
    
    tr_mat = np.dot(move1_matrix, rot_matrix)
    tr_mat = np.dot(tr_mat, flip_matrix)
    tr_mat = np.dot(tr_mat, move2_matrix)
    #print(tr_mat)
    matrix_elements = [tr_mat[0,0], tr_mat[0,1], tr_mat[1,0], tr_mat[1,1], tr_mat[0,2], tr_mat[1,2]]
    return affinity.affine_transform(polygon, matrix_elements)

def rotate_flip_all_polygons(ROIs, image, rot_angle, flipX, flipY):
    New_polygons = []
    for roi in ROIs:
        polygon = rotate_flip_polygon(Polygon(roi['points']), [image.getSizeX()/2, image.getSizeY()/2], rot_angle, flipX, flipY)
        New_polygons.append(polygon)
    return New_polygons

def read_tissue_positions_SR(spaceranger_path):
    # spaceranger v1.30
    tissue_positions = os.path.join(spaceranger_path,"spatial/tissue_positions_list.csv")
    if not os.path.exists(tissue_positions):
        # spacernger v2.10
        tissue_positions = os.path.join(spaceranger_path,"spatial/tissue_positions.csv")
    
    log.info(f"Reading spots coordintes from '{tissue_positions}'.")
    try:
        # spacernger v2.10
        # spatial/tissue_positions.csv has a header row
        df = pd.read_csv(
            tissue_positions,
            index_col="barcode"
        )
    except ValueError:
        # spaceranger v1.30
        # spatial/tissue_positions_list.csv doesn't have a header row
        df = pd.read_csv(
            tissue_positions,
            index_col="barcode",
            names=["barcode","in_tissue","array_row","array_col","pxl_row_in_fullres","pxl_col_in_fullres"]
       )
        # let's make sure pixel coordinates are numbers and in_tisssue is boolean
    df = df.astype({"in_tissue": bool, "pxl_row_in_fullres": int, "pxl_col_in_fullres":int})
    # filter spots that are only in the tissue
    df_in_tissue = df[df['in_tissue']]
    log.info(f"Total spots: {len(df)}. Spots in_tissue: {len(df_in_tissue)}.")
    return df_in_tissue

def assign_barcode_to_annotation(df_in_tisssue, ROIs, New_polygons):
    annotations = []
    log.info(f"Calling Spots inside ROIs.")
    for barcode,row in tqdm(df_in_tisssue.iterrows(), total=df_in_tisssue.shape[0]):
        roi_annoation = []
        #for roi in ROIs:
        for roi, polygon in zip(ROIs, New_polygons):
            #print(roi)
            #polygon = Polygon(roi['points'])
            point = Point([row['pxl_col_in_fullres'],row['pxl_row_in_fullres']])
            if polygon.contains(point):
                roi_annoation.append(roi['name'])
        if not roi_annoation:
            roi_annoation = ["N/A"]
        # why ";".join? because some annotations may overlap for the same spot/barcode
        # in that case we capture all and merge them with ';' 
        annotations.append([ barcode, ";".join(set(roi_annoation)) ])
    df_annotations = pd.DataFrame(annotations, columns=["barcode","ROI"]).set_index("barcode")

    return df_annotations, annotations

def plot_small_image(adata, folder_out, sample_name):
    sample_path = folder_out + '/' + sample_name
    plots = sq.pl.spatial_scatter(adata,color="ROI",alpha=0.75,size=1.5,save = str(sample_name))

def save_barcodes_ann_csv(adata, folder_out, sample_name):
    final_table = adata.obs
    final_table = final_table.drop(['in_tissue'], axis = 1)
    filename = folder_out + '/' + str(sample_name) + '_barcodes_ROIs.csv'
    final_table.to_csv(filename)

def main(csv_path, out_folder):
    table_input = pd.read_csv(csv_path)
    #initialization
    omero_host, omero_username, omero_password = get_OMERO_credentials()

    for i in range(table_input.shape[0]):
        
        spaceranger_path = table_input['Path'][i]
        omero_image_id = table_input['ImageID'][i]
        sample_name = table_input['Sample'][i]
        rot_angle = float(table_input['Rotation'][i])
        flipX = bool(table_input['FlipHorizontal'][i])
        flipY = bool(table_input['FlipVertical'][i])
        print('Sample ID: ' + str(omero_image_id))
        ROIs, image = collect_ROIs_from_OMERO(omero_username, omero_password, omero_host, omero_image_id)
        adata = read_SR_to_anndata(spaceranger_path)
        df_in_tissue = read_tissue_positions_SR(spaceranger_path)
        
        
        New_polygons = rotate_flip_all_polygons(ROIs, image, rot_angle, flipX, flipY)
        df_annotations, annotations = assign_barcode_to_annotation(df_in_tissue, ROIs, New_polygons)
        assert len(adata.obs)==len(df_annotations), "Inconsistent length between observations and annotations"
        adata.obs = pd.concat([adata.obs, df_annotations],axis=1)
        adata.obs['ROI'] = adata.obs['ROI'].astype('category')
        
        plot_small_image(adata, out_folder, sample_name)
        save_barcodes_ann_csv(adata, out_folder, sample_name)

if __name__ == "__main__":
    fire.Fire(main) 
