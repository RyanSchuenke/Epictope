from argparse import ArgumentParser
from epictope.plot_scores import plot_scores

def add_parser(subparsers:ArgumentParser):
    parser = subparsers.add_parser("plot", help="Plot the moving average 'min' score from a score file output produced by `epictope score`")
    
    parser.add_argument(dest="score_file", help="Score file csv produced by `epictope score`", type=str)

    parser.add_argument("-c", "--components", dest="components", action="store_true", required = False, 
                    help = "plot the moving averages of the individual score components rather than the default minimum score",
                    default = False,)

    parser.add_argument("-o", "--output", dest="output_file", required = False, 
                    help = "output file path for plotted score",
                    default = "",)
    
    parser.add_argument("-w", "--window-size", dest="window_size", required = False, 
                        help = "size of the sliding window used to average the EpicTope min score",
                        type = int, default = 7,)
    
    parser.add_argument("-n", "--non-zero", dest="min_non_zero_count", required = False, 
                        help = "the number of non-zero min scores needed for the average score to be calculated, othwise set to zero",
                        type = int, default = 1, )
    
    parser.set_defaults(func=run)

def run(args) -> None:
    if args.components:
        plot_scores(scores_file=args.score_file, output_file=args.output_file, 
                    window_size=args.window_size, min_non_zero_count=args.min_non_zero_count, 
                    score_components=["rsa", "ss_score", "inv_anchor2", "normalized_entropy"])
    else:
        plot_scores(scores_file=args.score_file, output_file=args.output_file, 
                    window_size=args.window_size, min_non_zero_count=args.min_non_zero_count)
