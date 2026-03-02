import argparse
from epictope.single_score import single_score
def add_parser(subparsers:argparse.ArgumentParser):
    parser = subparsers.add_parser("score", help="Calculate the scores for epictope site insertion for a protein from Uniprot", 
                                   formatter_class=argparse.RawTextHelpFormatter, 
                                   epilog="a local installation of iupred2a can be accessed in the working directory or provided with the environment variable 'IUPRED2A_PATH'.")
    
    parser.add_argument(dest="query", help="Uniprot accession or fasta file path of query sequence. If fasta file used, a custom structure is also required.", type=str)

    parser.add_argument("-p", "--plot", dest="plot", action="store_true", required = False, 
                    help = "flag for whether the outputted score file should be plotted"
                    )

    parser.add_argument("-c", "--config", dest="config", required = False, 
                    help = """path to config.yml file with parameter values for:
  'species': list of species to use as blast databases more msa and shannon entropy calculation
  'weights': dictionary of score component weights for shannon entropy 'h_weight', relative solvent accessibility 'rsa_weight', secondary structure 'ss_weight', and disordered binding region 'br_weight'
  'ss_key': dictionary with the secondary structure mapping for score of secondary structure characters [GHIECTBSP-] from dssp
  'max_sasa': dictionary mapping the maximum solvent accessibility for each amino acid""",
                    default = None)

    custom_struct_group = parser.add_argument_group(title='custom_structure', description='Arguments for custom structure input')

    custom_struct_group.add_argument("-s", "--structure", dest="custom_struct", required=False,
                    help="custom cif or pdb structure file. Currently tolerates missing the beginning or end of the protein, but not having gaps in the middle of the structure",
                    default = None,)
    
    parser.set_defaults(func=run)

def run(args) -> None:
    single_score(query=args.query, config_path=args.config, custom_struct=args.custom_struct, plot=args.plot)
