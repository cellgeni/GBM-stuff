import os
import io
import json
import logging
from getpass import getpass
from numpy import random

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
import yaml

log = logging.getLogger(__name__)


def get_OMERO_credentials():
    logging.basicConfig(format='%(asctime)s %(message)s')    
    log.setLevel(logging.DEBUG)
    omero_host = "wsi-omero-prod-02.internal.sanger.ac.uk"
    omero_username = input("Username$")
    omero_password = getpass("Password$")
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


def assign_cell_to_annotation(segmentation_table, ROIs, New_polygons, column_name_x, column_name_y, pixelsize):
    annotations = []
    log.info(f"Calling cell positions inside ROIs.")
    for barcode,row in tqdm(segmentation_table.iterrows(), total=segmentation_table.shape[0]):
        roi_annoation = []
        #for roi in ROIs:
        for roi, polygon in zip(ROIs, New_polygons):
            #print(roi)
            #polygon = Polygon(roi['points'])
            point = Point([row[column_name_x]/pixelsize,row[column_name_y]/pixelsize])
            if polygon.contains(point):
                roi_annoation.append(roi['name'])
        if not roi_annoation:
            roi_annoation = ["N-A"]
        # why ";".join? because some annotations may overlap for the same spot/barcode
        # in that case we capture all and merge them with ';' 
        annotations.append([ barcode, "; ".join(set(roi_annoation)) ])
    return annotations


def ReadConfFile(FilePath):
    with open(FilePath, 'r') as file:
        data = yaml.safe_load(file)
    omero_image_id = data['omero_image_id']
    column_name_x = data['column_name_x']
    column_name_y = data['column_name_y']
    rot_angle = data['rot_angle']
    flipX = data['flipX']; flipY = data['flipY']
    pixelsize = data['pixelsize']
    out_folder = data['output_folder']
    save_images = data['save_images']
    segmentation_csv_path_list = []
    for p in data['segmentation_csv']:
        segmentation_csv_path_list.append(p)
    return omero_image_id, column_name_x, column_name_y, segmentation_csv_path_list, rot_angle, flipX, flipY, pixelsize, out_folder, save_images

def plot_small_image(segmentation_table, out_folder, sample_name, column_name_x, column_name_y):
    i=0
    ann_list_unique = list(set(segmentation_table['annotation']))
    fig, axs = plt.subplots(1, 1, figsize=(20, 20))
    for name_roi in ann_list_unique:
        #index_list = list(filter(lambda x: ann_list[x] == name_roi, range(len(ann_list))))
        sub_segm_table = segmentation_table[segmentation_table['annotation']==name_roi]
        random_color =  [random.randint(100)/100, random.randint(100)/100, random.randint(100)/100]
        plt.scatter(sub_segm_table[column_name_x], sub_segm_table[column_name_y], color = random_color, s =1)
    path_fig = out_folder + '/' + sample_name[:-2] + '.png'
    plt.axis('equal')
    plt.gca().invert_yaxis()
    axs.legend(ann_list_unique, fontsize = 15, markerscale = 5)
    plt.savefig(path_fig, format="jpg", dpi = 300)


def add_annotations_to_table(segmentation_table, annotations):
    ann_list = []
    for a in annotations:
        ann_list.append(a[1])
    segmentation_table['annotation'] = ann_list
    
    return segmentation_table


def main(ConfFilePath):
    
    #initialization
    omero_image_id, column_name_x, column_name_y, segmentation_csv_path_list, rot_angle, flipX, flipY, pixelsize, out_folder, save_images = ReadConfFile(ConfFilePath)
    omero_host, omero_username, omero_password = get_OMERO_credentials()
    ROIs, image = collect_ROIs_from_OMERO(omero_username, omero_password, omero_host, omero_image_id)

    if len(ROIs) == 0:
        print('0 ROIs were found! Check whether you have any annotations in OMERO or are you owner of dataset?')
    else:
        for segmentation_csv_path in segmentation_csv_path_list:
            sample_name = os.path.basename(segmentation_csv_path)
            print('Working on segnmentation for: ' + sample_name)
            segmentation_table = pd.read_csv(segmentation_csv_path)
            if rot_angle!=0 and flipX!=False and flipY!=False:
                New_polygons = rotate_flip_all_polygons(ROIs, image, rot_angle, flipX, flipY)
            else:
                New_polygons = []
                for roi in ROIs:
                    New_polygons.append(Polygon(roi['points']))
            annotations = assign_cell_to_annotation(segmentation_table, ROIs, New_polygons, column_name_x, column_name_y, pixelsize)
            segmentation_table = add_annotations_to_table(segmentation_table, annotations)
            path_csv = out_folder + '/' +  sample_name
            segmentation_table.to_csv(path_csv)
            if save_images:
                plot_small_image(segmentation_table, out_folder, sample_name, column_name_x, column_name_y)
        
        

if __name__ == "__main__":
    fire.Fire(main) 
