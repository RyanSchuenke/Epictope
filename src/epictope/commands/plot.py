from argparse import ArgumentParser
from epictope.plot_scores import plot_scores, plot_score_components

def add_parser(subparsers:ArgumentParser):
    parser = subparsers.add_parser("plot", help="Plot the moving average 'min' score from a score file output produced by `epictope score`")
    
    parser.add_argument(dest="score_file", help="Score file csv produced by `epictope score`", type=str)

    parser.add_argument("-c", "--components", dest="components", action="store_true", required = False, 
                    help = "plot the moving averages of the individual score components rather than the default minimum score",
                    default = False,)

    parser.add_argument("-o", "--output", dest="output_file", required = False, 
                    help = "output file path for plotted score",
                    default = "",)
    
    parser.set_defaults(func=run)

def run(args) -> None:
    if args.components:
        plot_score_components(scores_file=args.score_file, output_file=args.output_file)
    else:
        plot_scores(scores_file=args.score_file, output_file=args.output_file)
