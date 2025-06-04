import os
import yaml

projname = 'fw-pGENMi'
SRCPATH = os.path.dirname(os.path.realpath(__file__))
PROJPATH = os.path.dirname(SRCPATH)

configs_dir = f'{PROJPATH}/configs'
results_dir = f'{PROJPATH}/results'
input_dir = f'{PROJPATH}/input'

def make_abspath(curr_path, append_path):
    # check if relative or abs path is provided
    # append_path makes the curr_path absolute
    if not os.path.isabs(curr_path):
        curr_path = os.path.realpath(os.path.join(append_path, curr_path))
    return curr_path

def read_yaml(file_path):
    with open(file_path) as stream:
        try:
            yaml_data = yaml.safe_load(stream)
            return yaml_data
        except yaml.YAMLError as exc:
            print(f"Error reading YAML file: {exc}")
            return None