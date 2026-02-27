import argparse
import sys
from epictope.commands.score import add_parser as add_score
from epictope.commands.plot import add_parser as add_plot
import logging
from datetime import datetime
import os

logger = logging.getLogger(__name__)

class EpictopeParser(argparse.ArgumentParser):
    def error(self, message):
        logger.error(message)
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
    os.makedirs("logs", exist_ok=True)
    
    # Create timestamp string
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

    # Create unique log filename
    if args.custom_struct:
        log_filename='logs/'+os.path.basename(args.query)+'_'+os.path.basename(args.custom_struct)+'_'+timestamp+'.log'
    else:
        log_filename='logs/'+os.path.basename(args.query)+'_'+timestamp+'.log'

    # Configure logging
    logging.basicConfig(
        filename=log_filename,
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(filename)s - %(message)s"
    )
    logger.info('Started')
    args.func(args)
    
    logger.info('Finished')
