''' The functions in this script have been adpated from, 
    https://knoweng.github.io 
    https://github.com/KnowEnG
'''

import os
import csv
import shutil
#from tempfile import mkdtemp
import urllib.request
from IPython.display import HTML

from kndatacleanup import data_cleanup
from knsamplesclustering import samples_clustering
from kngeneralclustering import general_clustering

from math import log10, floor

#from bs4 import BeautifulSoup
from html.parser import HTMLParser

from fwpgenmi.io_cfg import results_dir, make_abspath

REDIS_PARAMS = {
    'host': 'knowredis.knoweng.org',
    'password': 'KnowEnG',
    'port': 6379
}

class HTMLFilter(HTMLParser):
    text = ""
    def handle_data(self, data):
        self.text += "\n"
        self.text += data

def round_to_1(x):
    return round(x, -int(floor(log10(abs(x)))))


def fetch_network(edge_file_path, NETWORK_DIR_PATH):
    """Given the local path to an edge file, ensures that the edge file exists
    on disk, downloading it from AWS if necessary.

    Arguments:
        edge_file_path (str): The local path to an edge file as found in the
            data returned by `get_interaction_networks` and
            `get_gene_property_networks`.

    Returns:
        None

    """
    if not os.path.isfile(edge_file_path):
        url = "https://s3.amazonaws.com/KnowNets/KN-20rep-1706/" + \
            "userKN-20rep-1706/" + edge_file_path[len(NETWORK_DIR_PATH):]
        os.makedirs(os.path.dirname(edge_file_path), exist_ok=True)
        with urllib.request.urlopen(url) as response:
            with open(edge_file_path, 'wb') as out_file:
                shutil.copyfileobj(response, out_file)

def fetch_network_metadata(NETWORK_DIR_PATH):
    """Downloads from AWS the network overview metadata files required to
    implement the network utility methods.

    Arguments:
        None

    Returns:
        None

    """
    filenames = ['db_contents.txt', 'species_desc.txt', 'edge_type.txt']
    for filename in filenames:
        out_file_path = os.path.join(NETWORK_DIR_PATH, filename)
        if not os.path.isfile(out_file_path):
            url = "https://s3.amazonaws.com/KnowNets/KN-20rep-1706/" + \
                "userKN-20rep-1706/" + filename
            with urllib.request.urlopen(url) as response:
                with open(out_file_path, 'wb') as out_file:
                    shutil.copyfileobj(response, out_file)

def get_path_to_newest_file_having_prefix(search_dir_path, prefix):
    """Finds all files in `search_dir_path` whose name begins with `prefix`.
    Returns the newest of these matching files, or returns None if no matching
    file exists.

    Arguments:
        search_dir_path (str): The local path to the directory to search.
        prefix (str): The string used to filter the files in `search_dir_path`.

    Returns:
        str: The path to the newest matching file, or None if no matching files
            exist.

    """
    matches = [os.path.join(search_dir_path, name) \
        for name in os.listdir(search_dir_path) \
        if name.startswith(prefix)]
    # ensure they're all files
    matches = [m for m in matches if os.path.isfile(m)]
    return_val = None
    if matches:
        return_val = sorted(matches, \
            key=lambda path: os.path.getctime(path), reverse=True)[0]
    return return_val

def get_cleaned_file_path(original_file_path, results_dir_path):
    """Given the name of a file passed to `kndatacleanup.data_cleanup`,
    along with the `results_dir_path` passed to `kndatacleanup.data_cleanup`,
    returns the path at which the cleaned version of the input can be found.

    Arguments:
        original_file_path (str): The path to the file that was passed to
            `kndatacleanup.data_cleanup`.
        results_dir_path (str): The path to the results directory that was
            passed to `kndatacleanup.data_cleanup`.

    Returns:
        str: The path to the cleaned version of `original_file_path` that was
            or would be produced by `kndatacleanup.data_cleanup`.

    """
    original_name = os.path.basename(original_file_path)
    original_name_root = os.path.splitext(original_name)[0]
    return os.path.join(results_dir_path, original_name_root + "_ETL.tsv")

def get_gene_map_file_path(original_file_path, results_dir_path):
    """Given the name of an omics file passed to `kndatacleanup.data_cleanup`,
    along with the `results_dir_path` passed to `kndatacleanup.data_cleanup`,
    returns the path at which the mapping of gene names to gene identifiers
    can be found.

    Arguments:
        original_file_path (str): The path to the omics file that was passed to
            `kndatacleanup.data_cleanup`.
        results_dir_path (str): The path to the results directory that was
            passed to `kndatacleanup.data_cleanup`.

    Returns:
        str: The path to the gene-name mapping file that was or would be
            produced by `kndatacleanup.data_cleanup`.

    """
    original_name = os.path.basename(original_file_path)
    original_name_root = os.path.splitext(original_name)[0]
    return os.path.join(results_dir_path, original_name_root + "_MAP.tsv")


def get_network_species(NETWORK_DIR_PATH):
    """Returns information about the species found in the knowledge network.

    Arguments:
        None

    Returns:
        list: A list in which each element is a dictionary. Each dictionary has
            keys 'id' (which can be passed to other methods that require a
            `species_id`), 'short_latin_name', 'latin_name', 'familiar_name',
            and 'group_name'.

    """
    return_val = []
    species_file_path = os.path.join(NETWORK_DIR_PATH, 'species_desc.txt')
    with open(species_file_path) as csvfile:
        for row in csv.reader(csvfile, delimiter='\t'):
            return_val.append({
                'id': row[0],
                'short_latin_name': row[1],
                'latin_name': row[2],
                'familiar_name': row[3],
                'group_name': row[5]
            })
    return return_val

def display_network_species():
    """Displays a table of the species found in the knowledge network.

    Arguments:
        None

    Returns:
        None

    """
    html_string = "<table><tr><th>Familiar Name (Latin Name)</th><th>Species Id</th></tr>"
    for species in get_network_species():
        html_string += "<tr><td>" + species['familiar_name'] + " (" + \
            species['latin_name'] + ")</td><td>" + species['id'] + "</td></tr>"
    html_string += "</table>"
    return HTML(html_string)

def get_edge_type_name_to_pretty_name(NETWORK_DIR_PATH):
    """Returns a dictionary in which the keys are edge type names and the values
    are pretty network names.

    Arguments:
        None

    Returns:
        dict: A dictionary in which the keys are edge type names and the values
            are pretty network names.

    """
    return_val = {}
    file_path = os.path.join(NETWORK_DIR_PATH, 'edge_type.txt')
    with open(file_path) as csvfile:
        for row in csv.DictReader(csvfile, delimiter='\t'):
            return_val[row['et_name']] = row['pretty_name']
    return return_val

def get_interaction_networks(species_id, NETWORK_DIR_PATH):
    """Given a `species_id`, returns information about the interaction networks
    available in the knowledge network.

    Arguments:
        species_id (int or str): The id for the species of interest, as returned
            by `get_network_species` or displayed by `display_network_species`.

    Returns:
        list: A list in which each element is a dictionary. Each dictionary has
            two keys, 'name' and 'edge_file_path'.

    """
    species_id = str(species_id) # user-friendliness
    return_val = []
    contents_file_path = os.path.join(NETWORK_DIR_PATH, 'db_contents.txt')
    with open(contents_file_path) as csvfile:
        edge_type_name_to_pretty_name = get_edge_type_name_to_pretty_name(NETWORK_DIR_PATH)
        for row in csv.DictReader(csvfile, delimiter='\t'):
            if row['n1_type'] == 'Gene' and row['taxon'] == species_id:
                return_val.append({
                    'name': edge_type_name_to_pretty_name[row['et_name']],
                    'edge_file_path': os.path.join(\
                        NETWORK_DIR_PATH, 'Gene', species_id, row['et_name'], \
                        species_id + '.' + row['et_name'] + '.edge')
                })
    return return_val

def display_interaction_networks(species_id, NETWORK_DIR_PATH):
    """Given a `species_id`, displays information about the interaction
    networks available in the knowledge network.

    Arguments:
        species_id (int or str): The id for the species of interest, as returned
            by `get_network_species` or displayed by `display_network_species`.

    Returns:
        None

    """
    html_string = "<table><tr><th>Interaction Network Name</th><th>Edge File Path</th></tr>"
    for network in get_interaction_networks(species_id, NETWORK_DIR_PATH):
        html_string += "<tr><td>" + network['name'] + "</td><td>" + \
            network['edge_file_path'] + "</td></tr>"
    html_string += "</table>"
    return HTML(html_string)

def get_gene_property_networks(species_id, NETWORK_DIR_PATH):
    """Given a `species_id`, returns information about the gene-property
    networks available in the knowledge network.

    Arguments:
        species_id (int or str): The id for the species of interest, as returned
            by `get_network_species` or displayed by `display_network_species`.

    Returns:
        list: A list in which each element is a dictionary. Each dictionary has
            two keys, 'name' and 'edge_file_path'.

    """
    species_id = str(species_id) # user-friendliness
    return_val = []
    contents_file_path = os.path.join(NETWORK_DIR_PATH, 'db_contents.txt')
    with open(contents_file_path) as csvfile:
        edge_type_name_to_pretty_name = get_edge_type_name_to_pretty_name(NETWORK_DIR_PATH)
        for row in csv.DictReader(csvfile, delimiter='\t'):
            if row['n1_type'] == 'Property' and row['taxon'] == species_id:
                return_val.append({
                    'name': edge_type_name_to_pretty_name[row['et_name']],
                    'edge_file_path': os.path.join(\
                        NETWORK_DIR_PATH, 'Property', species_id, row['et_name'], \
                        species_id + '.' + row['et_name'] + '.edge')
                })
    return return_val

def display_gene_property_networks(species_id, NETWORK_DIR_PATH):
    """Given a `species_id`, displays information about the gene-property
    networks available in the knowledge network.

    Arguments:
        species_id (int or str): The id for the species of interest, as returned
            by `get_network_species` or displayed by `display_network_species`.

    Returns:
        None

    """
    html_string = "<table><tr><th>Interaction Network Name</th><th>Edge File Path</th></tr>"
    for network in get_gene_property_networks(species_id, NETWORK_DIR_PATH):
        html_string += "<tr><td>" + network['name'] + "</td><td>" + \
            network['edge_file_path'] + "</td></tr>"
    html_string += "</table>"
    return HTML(html_string), html_string


def do_clustering(\
    omics_file_path, phenotype_file_path, results_dir_path, num_clusters, \
    species_id, interaction_network_edge_file_path, network_influence, \
    num_bootstraps, bootstrap_sample_fraction, NUM_CPUS= 2):
    """Performs a clustering upon the samples found in `omics_file_path`.

    Arguments:
        omics_file_path (str): The path to the omics file.
        phenotype_file_path (str): The path to a file containing phenotype data
            on the same samples as found in `omics_file_path`, or None if no
            phenotype data are to be analyzed. If analyzed, each phenotype will
            scored for statistically significant differences between the
            clusters.
        results_dir_path (str): The path to a directory where results files
            should be stored.
        num_clusters (int): The number of clusters to create.
        species_id (int or str): The id for the species of interest, as returned
            by `get_network_species` or displayed by `display_network_species`,
            or None not using an `interaction_network_edge_file_path`.
        interaction_network_edge_file_path (str): The path to an interaction
            network edge file, to use a knowledge-guided approach to clustering,
            or else None.
        network_influence (float): A number between 0 and 1 that specifies the
            amount to which network data should influence the results, or None
            if not using an `interaction_network_edge_file_path`.
        num_bootstraps (int): A number of bootstrap iterations to run. Use 0 for
            no bootstrapping.
        bootstrap_sample_fraction (float): A number between 0 and 1 that
            specifies what fraction of the data should be used in each bootstrap
            iteration, or None if not using bootstrapping.

    Returns:
        None

    """
    try:
        species_id = str(species_id) # user-friendliness
        os.makedirs(results_dir_path, exist_ok=True)

        if interaction_network_edge_file_path is None:
            pipeline_type = 'general_clustering_pipeline'
        else:
            #fetch_network(interaction_network_edge_file_path)
            pipeline_type = 'samples_clustering_pipeline'

        cleanup_parameters = {
            'spreadsheet_name_full_path': omics_file_path,
            'pipeline_type': pipeline_type,
        'results_directory': results_dir_path
        }
        if phenotype_file_path is not None:
            cleanup_parameters['phenotype_name_full_path'] = phenotype_file_path
        if interaction_network_edge_file_path is not None:
            cleanup_parameters.update({
                'gg_network_name_full_path': interaction_network_edge_file_path,
                'taxonid': species_id,
                'source_hint': '',
                'redis_credential': {
                    'host': REDIS_PARAMS['host'],
                    'port': REDIS_PARAMS['port'],
                    'password': REDIS_PARAMS['password']
                }
            })
        data_cleanup.run_pipelines(cleanup_parameters, \
            data_cleanup.SELECT[pipeline_type])

        clustering_parameters = {
            'spreadsheet_name_full_path': get_cleaned_file_path(\
                omics_file_path, results_dir_path),
            'results_directory': results_dir_path,
            'processing_method': 'parallel',
            'parallelism': NUM_CPUS,
            'number_of_clusters': num_clusters,
            'run_directory': results_dir_path,
            'tmp_directory': './tmp'
        }
        if phenotype_file_path is not None:
            clustering_parameters.update({
                'phenotype_name_full_path': get_cleaned_file_path(\
                    phenotype_file_path, results_dir_path),
                'threshold': 15
            })

        method_prefix = ''
        if num_bootstraps > 0:
            clustering_parameters.update({
                'number_of_bootstraps': num_bootstraps,
                'rows_sampling_fraction': 1.0,
                'cols_sampling_fraction': bootstrap_sample_fraction
            })
            method_prefix = 'cc_'

        if interaction_network_edge_file_path is not None:
            clustering_parameters.update({
                'gg_network_name_full_path': interaction_network_edge_file_path,
                'rwr_max_iterations': 100,
                'rwr_convergence_tolerence': 1.0e-4,
                'rwr_restart_probability': network_influence,
                'top_number_of_genes': 100,
                'nmf_conv_check_freq': 50,
                'nmf_max_invariance': 200,
                'nmf_max_iterations': 10000,
                'nmf_penalty_parameter': 1400,
                'method': method_prefix + 'net_nmf'
            })
            samples_clustering.SELECT[clustering_parameters['method']](\
                clustering_parameters)
        else:
            clustering_parameters.update({
                'top_number_of_rows': 100,
                'affinity_metric': 'euclidean',
                'linkage_criterion': 'ward',
                'method': method_prefix + 'hclust'
            })
            general_clustering.SELECT[clustering_parameters['method']](\
                clustering_parameters)
    except:
        print("Something went wrong! Check the debugging information below, " + \
            "and look for log output in " + results_dir_path)
        raise
    else:
        print("Find results in " + results_dir_path)
        

        
def main(NETWORK_DIR_PATH, species_id= 9606):
    # fetch the network metadata and get the networks available for each species 
    # you can find the species id from the network meta data
    os.makedirs(NETWORK_DIR_PATH, exist_ok=True)
    fetch_network_metadata(NETWORK_DIR_PATH)
    _, html_string= display_gene_property_networks(species_id, NETWORK_DIR_PATH)
    f = HTMLFilter()
    f.feed(html_string)
    file= open(f"{NETWORK_DIR_PATH}/gene_nw_{species_id}", "w")
    file.write(f.text)
    file.close()
    
    #for interaction_network_edge_file_path in ['COCA_results/network/Property/9606/PPI_complex/9606.PPI_complex.edge']:
    #    fetch_network(interaction_network_edge_file_path, NETWORK_DIR_PATH)
    #print(f.text)

    
if __name__ == "__main__":
    # execute only if run as a script
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--NETWORK_DIR_PATH', default= './network/', type=str)
    parser.add_argument('--species_id', default= 9606, type=int)
    args = parser.parse_args()
    
    # check for absolute vs relative path
    NETWORK_DIR_PATH = args.NETWORK_DIR_PATH
    NETWORK_DIR_PATH = make_abspath(NETWORK_DIR_PATH, results_dir)
    main(NETWORK_DIR_PATH, args.species_id)
    
