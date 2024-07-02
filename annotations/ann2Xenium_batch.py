import fire
import yaml
import ann2Xenium_v2
import logging
log = logging.getLogger(__name__)

def ReadConfFile(FilePath):
    with open(FilePath, 'r') as file:
        data = yaml.safe_load(file)
    
    
    omero_username = data['omero_username']
    omero_password = data['omero_password']
    
    conf_file_path_list = []
    for p in data['conf_file_paths']:
        conf_file_path_list.append(p)
    
    return conf_file_path_list, omero_username, omero_password


def get_OMERO_credentials():
    logging.basicConfig(format='%(asctime)s %(message)s')    
    log.setLevel(logging.DEBUG)
    
    omero_username = input("Username$")
    omero_password = getpass("Password$")
    return omero_username, omero_password


def main(ConfFilePathMaster):
    
    omero_host = "wsi-omero-prod-02.internal.sanger.ac.uk"
    conf_file_path_list, omero_username, omero_password = ReadConfFile(ConfFilePathMaster)
    
    if omero_username == 'None' or omero_password == 'None':
        omero_username, omero_password = get_OMERO_credentials()
    
   
    
    for conf_file_path in conf_file_path_list:
        ann2Xenium_v2.main(conf_file_path, omero_username, omero_password, omero_host)
    
    
    
if __name__ == "__main__":
    fire.Fire(main) 


