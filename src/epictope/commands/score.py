from argparse import ArgumentParser
from epictope.single_score import single_score
def add_parser(subparsers:ArgumentParser):
    parser = subparsers.add_parser("score", help="Calculate the scores for epictope site insertion for a protein from Uniprot")
    
    parser.add_argument(dest="query", help="Uniprot accession of query sequence", type=str)
    parser.add_argument("-p", "--plot", dest="plot", action="store_true", required = False, 
                    help = "boolean value for if the "
                    )

    parser.add_argument("-c", "--config", dest="config", required = False, 
                    help = "path to config.yml file with parameter values",
                    default = None)

    custom_struct_group = parser.add_argument_group(title='custom_structure', description='Arguments for custom structure input')

    custom_struct_group.add_argument("-s", "--structure", dest="custom_struct", required=False,
                    help="custom cif or pdb structure file",
                    default = None,)

    custom_struct_group.add_argument("-i", "--index", dest="start_index", required=False,
                    help="starting amino acid position in structure file, 1 indexed",
                    default=1, type=int)
    
    parser.set_defaults(func=run)

def run(args) -> None:
    single_score(query=args.query.upper(), config_path=args.config, custom_cif=args.custom_struct, res_start=args.start_index, graph=args.plot)
