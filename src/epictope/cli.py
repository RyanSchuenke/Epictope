import argparse
import sys
from epictope.commands.score import add_parser as add_score
from epictope.commands.plot import add_parser as add_plot

class EpictopeParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

def main() -> None:
    """
    Command line interface for running epictope
    """
    parser = EpictopeParser(prog="epictope",
        description= 'Runs the Epictope pipeline on the provided protein accession', 
        formatter_class=argparse.RawTextHelpFormatter)
    
    subparsers = parser.add_subparsers(
        title="commands",
        dest="command",
        required=True
    )
    add_score(subparsers)
    add_plot(subparsers)
    
    args = parser.parse_args()
    args.func(args)
