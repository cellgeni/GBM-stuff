import os
import io
import json
import logging
from getpass import getpass
import numpy as np
# conda create -n GBM python=3.10 -c conda-forge
# conda install scanpy squidpy pandas pillow shapely omero-py -c conda-forge
import omero
import scanpy as sc
import squidpy as sq
import pandas as pd
from tqdm import tqdm
from PIL import Image
from omero.gateway import BlitzGateway
from shapely import affinity,plotting
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
import os
import squidpy as sq
import fire


log = logging.getLogger(__name__)

def get_OMERO_credentials():
    logging.basicConfig(format='%(asctime)s %(message)s')    
    log.setLevel(logging.DEBUG)
    omero_host = "wsi-omero-prod-02.internal.sanger.ac.uk"
    omero_username = input("Username$")
    omero_password = getpass("Password$")
    return omero_host, omero_username, omero_password

def collect_ROIs_from_OMERO(omero_username, omero_password, omero_host, omero_image_id, path_ann_csv):
    ROIs = []
    log.info(f"Connecting to OMERO at {omero_host}")
    with BlitzGateway(omero_username, omero_password, host=omero_host, port=4064, secure=True) as conn:
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

            primary_shape = roi.getPrimaryShape()
            name = primary_shape.getTextValue().val
            #if ROI name is empty - renamed to "Non_labelled", makes sure first letter is capital, if additional csv provided - unified all different ROIs name into one
            name = rename_ROI(name, path_ann_csv)
            #print(name)
            #, separates x and y and  space separates points
            try:
                points = [(lambda xy : list(map(float,xy.split(","))))(xy) for xy in primary_shape.getPoints().val.split(" ")]
                ROIs.append({
                "name": name,
                "points": points
                })
            except:
                if primary_shape.__class__.__name__ == 'RectangleI':
                    points = get_corners_rectangle(primary_shape)
                    ROIs.append({
                    "name": name,
                    "points": points
                    })
                else:
                    pass
            
            log.debug(f"Found ROI id={roi.id.val} name='{name}' type={primary_shape.__class__.__name__}")

    log.info(f"Found {len(ROIs)} ROIs in total")
    return ROIs, image

def get_corners_rectangle(primary_shape):
    x0 = primary_shape.getX()._val; y0 = primary_shape.getY()._val
    w = primary_shape.getWidth()._val; h = primary_shape.getHeight()._val
    return [(x0, y0), (x0+w,y0), (x0+w,y0+h), (x0,y0+h)]

def rename_ROI(old_name, path_ann_name_csv=None):
    name = old_name
    if name == '': name = 'Not_Labelled'
    name = name.replace('/', '-')
    name = make_first_letter_upper(name)
    if path_ann_name_csv:
        ann_name = pd.read_csv(path_ann_name_csv)
        if (ann_name[ann_name['original'].isin([name])]).shape[0]>0:
            index = ann_name[ann_name['original'].isin([name])].index[0]
            name = ann_name['unified'][index]
    return str(name)
        
        

    
def make_first_letter_upper(some_string):
    if some_string[0].isupper() == False:
        some_string = some_string[0].upper() + some_string[1:]
    return some_string
 
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

def group_ROIs_by_group(ROIs, New_polygons):
    gROIs = {}
    for roi,pol in zip(ROIs, New_polygons):
        polygon = pol
        key = roi['name']
        if key in gROIs:
            gROIs[key].append(polygon)
        else:
            gROIs[key] = [polygon]
    return gROIs


def assign_barcode_rois(gROIs, df_in_tisssue, spaceranger_path):
    #read spot information
    with open(os.path.join(spaceranger_path,"spatial","scalefactors_json.json"),"rt") as f:
        scalefactors = json.load(f)
    spot_radius = scalefactors['spot_diameter_fullres']/2.0
    
    annotations = []
    #log.info(f"Calling {len(df_in_tisssue)} Spots inside {len(ROIs)} ROIs.")
    for barcode,row in tqdm(df_in_tisssue.iterrows(), total=df_in_tisssue.shape[0]):
        spot_annotations = {}
        
        spot = Point([row['pxl_col_in_fullres'],row['pxl_row_in_fullres']]).buffer(spot_radius)
        for roi_name,values in gROIs.items():
            coverage = 0.0
            for roi in values:
                try:
                    # to fix TopologyException: Self-intersection
                    roi = roi.buffer(0)
                    # if they intersect
                    if spot.intersects(roi):
                        # estimate how much they do intersect in percentage of spot area
                        iarea = spot.intersection(roi).area/spot.area
                        if iarea>coverage:
                            coverage = iarea
                except Exception as e:
                    log.error(f"Spot intersection failed. barcode={barcode} roi='{roi_name}'. Error: {repr(e)}")
            spot_annotations[roi_name] = coverage
        annotations.append({"barcode": barcode, "annotations":spot_annotations})
    data = [{"barcode": a["barcode"], **a["annotations"]} for a in annotations]
    df = pd.DataFrame(data).set_index("barcode")
    return df, spot_radius

def get_dict_with_parents_child(gROIs):
    hierarchy={}
    for pname,ppolygon in gROIs.items():
        hierarchy[pname] = {"children": [],"parents": []}
        for p in ppolygon:
            #print('New polygon of ' + pname)
            # check against all the other catergoires/roi names
            for cname,cpolygon in gROIs.items():
                # make sure it's not the same as the one we are checking for
                if cname == pname: continue
                for c in cpolygon:
                    # child test
                    if c.within(p):
                        #check if not in the list and add
                        if cname not in hierarchy[pname]["children"]:
                            hierarchy[pname]["children"].append(cname)
                    # paternity test
                    elif c.contains_properly(p):
                        if cname not in hierarchy[pname]["parents"]:
                            hierarchy[pname]["parents"].append(cname)

    return hierarchy


def get_changed_rois(dictionary):
    changed_rois = []
    for n_roi in dictionary.keys(): 
        if dictionary[n_roi]['changed']>0:
            changed_rois.append(n_roi)
    return changed_rois

def get_area_of_corresponding_polygon(polygon_list, point_x, point_y, spot_radius):
    spot = Point([point_x,point_y]).buffer(spot_radius)
    area_roi = np.nan
    for polygon in polygon_list:
        if spot.intersects(polygon):
            area_roi = polygon.area
            break
    return area_roi



def get_rois_non_zero(df_row):
    qq = np.where(df_row>0)
    list_of_rois = []
    if len(qq[1])>0:
        for i in range(len(qq[1])):
            list_of_rois.append(df_row.columns[qq[1][i]])
    else:
        list_of_rois = [None]
    return list_of_rois

def define_one_ROI_per_spot(df_rois, rois_dict, gROIs, spot_radius, df_in_tissue):
    print(df_rois.shape)
    print(df_in_tissue.shape)
    list_of_indeces = df_rois.index
    n = len(list_of_indeces)
    lst_Nones = [None for _ in range(n)]
    single_column_df = pd.DataFrame(index = list_of_indeces)
    single_column_df['ROI_one'] = lst_Nones

    
    single_column_df['ROI_one'] = single_column_df['ROI_one'].astype("string")
    for i in range(df_rois.shape[0]):
        row = df_rois.iloc[[i]]
        list_of_rois = get_rois_non_zero(row)
        #roi_name_max_level = 'Na
        if len(list_of_rois) == 1:
            
            if list_of_rois[0] is None:
                single_column_df['ROI_one'].iloc[i] = 'nan'
            else:
                single_column_df['ROI_one'].iloc[i] = list_of_rois[0]
        
        else:
            #search for the roi with the highest level
            if np.isnan(rois_dict[list_of_rois[0]]['level']):
                # if levels of annoations ROIs are not defined, I chose the smallest one by area 
                for roi_name in list_of_rois:
                    polygon_list = gROIs[roi_name]
                    pos_spot_x = df_in_tissue['pxl_col_in_fullres'][i]; pos_spot_y = df_in_tissue['pxl_row_in_fullres'][i];
                    area_polygon = get_area_of_corresponding_polygon(polygon_list, pos_spot_x, pos_spot_y, spot_radius)
                    print(roi_name)
                    print(area_polygon)
                    if roi_name==list_of_rois[0]:
                        area_min = area_polygon
                        roi_name_max_level = roi_name
                    else:
                        if area_polygon<area_min:
                            area_min = area_polygon
                            roi_name_max_level = roi_name
                    
                
            else:
                level_max = -1
                for roi_name in list_of_rois:
                    #print(roi_name)
                    if rois_dict[roi_name]['level']>level_max:
                        level_max = rois_dict[roi_name]['level']
                        roi_name_max_level = roi_name
            
            single_column_df['ROI_one'].iloc[i] = roi_name_max_level
            
    return single_column_df

def add_nlevel(dictionary):
    for n_roi in dictionary.keys(): 
        dictionary[n_roi]['level']=0
        dictionary[n_roi]['changed']=0

    for n_roi in dictionary.keys():
        parents = dictionary[n_roi]['parents']
        if len(parents)>0:
            n_lvl = []
            for parent_name in parents:
                n_lvl.append(dictionary[parent_name]['level'] + 1)
            dictionary[n_roi]['level'] = np.max(np.array(n_lvl))
            dictionary[n_roi]['changed'] = 1
    qq = get_changed_rois(dictionary)
    n_iter = 0
    while len(get_changed_rois(dictionary))>0 and n_iter<100:
        n_iter+=1
        for n_roi in dictionary.keys():
            parents = dictionary[n_roi]['parents']
            if len(parents)>0:
                n_lvl = []; 
                for parent_name in parents:
                    if dictionary[parent_name]['changed'] == 1:
                        n_lvl.append(dictionary[parent_name]['level'] + 1)
                if len(n_lvl)>0:
                    dictionary[n_roi]['level'] = np.max(np.array(n_lvl))
                    dictionary[n_roi]['changed'] = 1
                else:
                    dictionary[n_roi]['changed'] = 0
    
    #check if there was a logical eroor in ROIs hierarchy: Region A is inside Region B, and at the same time Region B is inside of region B
    if n_iter == 100:
        print('Impossible to determine ROis levels as there is a hierarchy error!')
        for n_roi in dictionary.keys():
            dictionary[n_roi]['level'] = np.nan

    for n_roi in dictionary.keys():
        del(dictionary[n_roi]['changed'])
    return dictionary

def plot_small_image(adata, folder_out, sample_name):
    sample_path = folder_out + '/' + sample_name
    try:
        plots = sq.pl.spatial_scatter(adata,color="ROI_one",alpha=0.75,size=2.5,save = str(sample_path))
    except:
        pass


def plot_images_rois(adata, folder_out, sample_name, dict_rois_level):
    for k in dict_rois_level.keys():
        image_path = folder_out + '/' + sample_name + '_' + k + '.png'
        try:
            sq.pl.spatial_scatter(adata, title = k, color=k, alpha=0.75, size=1.5, save = str(image_path))
        except:
            pass


def main(csv_path, out_folder, path_ann_csv = None, save_small_image = True, save_images_rois = True, save_csv = True):
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
        ROIs, image = collect_ROIs_from_OMERO(omero_username, omero_password, omero_host, omero_image_id, path_ann_csv)
        adata = read_SR_to_anndata(spaceranger_path)
        df_in_tissue = read_tissue_positions_SR(spaceranger_path)
        
        #rotate and flip polygons of ROIs
        New_polygons = rotate_flip_all_polygons(ROIs, image, rot_angle, flipX, flipY)

        #assign annotations and get gierarchy information about ROIs
        gROIs = group_ROIs_by_group(ROIs, New_polygons) # I update gROIs with also new polygons after rotation
        df_annotations, spot_radius = assign_barcode_rois(gROIs, df_in_tissue, spaceranger_path)
        #df_annotations = replace_symbol_by_dash_in_table(df_annotations)
        dict_rois = get_dict_with_parents_child(gROIs)
        
        dict_rois_level = add_nlevel(dict_rois)
        #dict_rois_level = replace_symbol_by_dash_in_dict(dict_rois_level)
        assert len(adata.obs)==len(df_annotations), "Inconsistent length between observations and annotations"
        
        #adding column with categorical annotation based on highest level roi
        single_col_df = define_one_ROI_per_spot(df_annotations, dict_rois_level, gROIs, spot_radius, df_in_tissue)
        df_annotations = pd.concat([df_annotations, single_col_df], axis = 1)
        if save_csv:
            fullpath = os.path.join(out_folder, sample_name + '.csv')
            df_annotations.to_csv(fullpath)
        
        #adding all information to anndata
        adata.obs = pd.concat([adata.obs, df_annotations],axis=1)
        
        adata.obs['ROI_one'] = adata.obs['ROI_one'].astype('object')
        adata.uns.update({'rois_hierarchy': dict_rois_level})
        print(adata.obs['ROI_one'])
        #save adata
        adata_path = out_folder + '/' + sample_name + '.h5ad'
        adata.write_h5ad(adata_path)
        
        #save images
        if save_small_image: 
            plot_small_image(adata, out_folder, sample_name)

        if save_images_rois: 
            plot_images_rois(adata, out_folder, sample_name, dict_rois_level)
        
        
        


if __name__ == "__main__":
    fire.Fire(main) 
