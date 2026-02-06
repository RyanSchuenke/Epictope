from argparse import ArgumentParser
from epictope.plot_scores import plot_scores

def add_parser(subparsers:ArgumentParser):
    parser = subparsers.add_parser("plot", help="Plot the score file output from `epictope score`")
    
    parser.add_argument(dest="score_file", help="Score file csv produced by `epictope score`", type=str)


    parser.add_argument("-o", "--output", dest="output_file", required = False, 
                    help = "output file path for plotted score",
                    default = "",)
    
    parser.set_defaults(func=run)

def run(args) -> None:
    plot_scores(scores_file=args.score_file, output_file=args.output_file)
