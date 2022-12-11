"""
### ======= Utility functions ======= ###
This snakemake file contains utility functions that are used throughout the pipeline. 

# Author: Will Hannon 
# Email: wwh22@uw.edu
"""

def aggregate_csv(csv_list, out_csv):
    """
    This function aggregates lists of CSVs from
    the input of a snakemake rule. It's intended
    to be used in any aggregation rules throughout 
    the pipeline. 

    Try/Except to deal with empty CSVs.
    """
    valid_csv_list = []
    for csv in csv_list:
        try:
            pd.read_csv(csv)
            valid_csv_list.append(csv)
        except: 
            pass
    df = pd.concat(map(pd.read_csv, valid_csv_list))
    df.to_csv(out_csv[0], index = False)



def get_readgroup(wildcards):
    """ 
    Get the sample name from the metadata file. This will be used as
    the readgroup for the alignment.
    """
    # Read the metadata into Pandas dataframe
    samples_df = pd.read_csv(config['samples']['file'])
    # Get the viral genome for a given accession
    readgroup = samples_df.loc[samples_df.Run == wildcards.accession, ["Sample"]].values.flatten().tolist()[0]
    # Return the readgroup
    return str(readgroup)