"""NeuN647 mouse LCT atlas"""

__version__ = "0"  # will be used to set minor version of the atlas

from bg_atlasgen.wrapup import wrapup_atlas_from_data
from bg_atlasgen.mesh_utils import *
import tifffile
from pathlib import Path
import pandas as pd
from tqdm import tqdm
import multiprocessing as mp
from functools import partial
from treelib import Tree
import os, shutil

def get_struct_path(df, structure_id_path):
    """
    Recurisive utility function for extracting the hierarchical structure path for a given brain region
    """
    
    if df.loc[df['id']==structure_id_path[-1]]['depth'].squeeze() == 0:
        structure_id_path.append(-1)
        structure_id_path.reverse()
        return structure_id_path
    else:
        structure_id_path.append(int(df.loc[df['id']==structure_id_path[-1]]['parent_structure_id'].squeeze()))
        return get_struct_path(df,structure_id_path)

def get_structures_list(structure_path, left_root_idx):
    """
    Create structure dictionary in line with what brainglobe expects
    """

    df = pd.read_csv(structure_path)
    df.loc[0,'id'] = left_root_idx # can't have id be 0 because it corresponds to black voxels
    df['parent_structure_id'] = df['parent_structure_id'].replace(0,left_root_idx)
    
    struct_list = [{"acronym":"root",
                    "name":"root",
                    "id":-1,
                    "structure_id_path":[-1],
                    "rgb_triplet":[255,255,255]}]
    for region_id in df['id']:
        df_ = df.loc[df['id']==region_id]
        structure_id_path = get_struct_path(df, [region_id])
        region_dict = {"acronym":df_['acronym'].squeeze(),
                        "name":df_['name'].squeeze(),
                        "id":region_id,
                        "structure_id_path": structure_id_path,
                        "rgb_triplet":[255,255,255]}
        struct_list.append(region_dict)
    return struct_list

def create_tree(struct_list):
    tree = Tree()
    tree.create_node('root',-1)
    for ele in struct_list:
        tree.create_node(ele['acronym'],ele['id'],parent=ele['structure_id_path'][-2])
    return tree

def create_atlas(working_dir, resolution, template_path, annotation_path, 
                 structure_path):
    """Function to generate source data for an atlas.

    Parameters
    ----------
    working_dir : Path object
        Path where atlas will be created.
    resolution :
        Resolution of the atlas, in um.

    Returns
    -------
    Path object
        Path to the final compressed atlas file.

    """

    ATLAS_NAME = "NeuN647_full"
    SPECIES = "Mus musculus"
    ATLAS_LINK = "https://lifecanvastech.com"
    CITATION = "unpublished"
    ORIENTATION = "sal" #sal = transverse orientation 
    
    if os.path.exists(working_dir):
        shutil.rmtree(working_dir)
        
    os.makedirs(working_dir)


    # do stuff to create the atlas
    # bg_atlasgen produces images 
    template_volume = tifffile.imread(template_path)  # volume with reference (read in as z,y,x)
    annotated_volume = tifffile.imread(annotation_path)  # volume with structures annotations
    structures_list = get_structures_list(structure_path, 20000)  # list of valid structure dictionaries
    root_id = -1 # id of the root structure
    
    # mesh creation -- takes too long and produces huge files
    # pth = Path(meshes_path)
    # tree = create_tree(structures_list)
    # labels = list(np.unique(annotated_volume))
    # meshes_dict = {} # dictionary of files with region meshes (struct_id : path)
    # print("Making region meshes...")
    # if num_workers:
    #     args_list = [(pth,node,tree,labels,annotated_volume,root_id,8,0.6,False) for node in tree.all_nodes()]
    #     p = mp.Pool(num_workers)
	# 	list(tqdm(p.imap(create_region_mesh,args_list), total = len(args_list)))
    #     p.close(); p.join()
    # else:
    #     for node in tqdm(tree.all_nodes()):
    #         args = (pth,node,tree,labels,annotated_volume,root_id,8,0.6,False)
    #         create_region_mesh(args)
    #         meshes_dict[node.identifier] = pth / f"{node.identifier}.obj"
    
    meshes_dict = None
    
    
    

    output_filename = wrapup_atlas_from_data(
        atlas_name=ATLAS_NAME,
        atlas_minor_version=__version__,
        citation=CITATION,
        atlas_link=ATLAS_LINK,
        species=SPECIES,
        resolution=(resolution,) * 3,  # if isotropic - highly recommended
        orientation=ORIENTATION,
        root_id=root_id,
        reference_stack=template_volume,
        annotation_stack=annotated_volume,
        structures_list=structures_list,
        meshes_dict=meshes_dict,
        working_dir=working_dir,
        hemispheres_stack=None,
        cleanup_files=False,
        compress=True,
    )

    return output_filename


if __name__ == "__main__":
    resolution = 10  # some resolution, in microns
    template_path = '/mnt/beegfs/ML_Data/Reg_Databases/atlas/NeuN647_10_full_transverse.tiff'
    annotation_path = '/mnt/beegfs/ML_Data/Reg_Databases/atlas/edge_annotation_10_full_transverse_LR.tiff'
    structure_path = '/mnt/beegfs/ML_Data/Reg_Databases/atlas/full_brain_regions_LR.csv'
    
    # Generated atlas path:
    bg_root_dir = "/home/user/.brainglobe"#"/home/user/Documents/brainreg_tests/atlases"

    create_atlas(bg_root_dir, resolution, template_path, annotation_path, structure_path)
