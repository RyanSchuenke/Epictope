import argparse
def main() -> None:
    parser = argparse.ArgumentParser(prog="epictope",
        description= 'Runs the Epictope pipeline on the provided protein accession')
    
    subparsers = parser.add_subparsers(
        title="commands",
        dest="command",
        required=True
    )
    from .commands.score import add_parser as add_score
    from .commands.plot import add_parser as add_plot
    add_score(subparsers)
    add_plot(subparsers)
    
    args = parser.parse_args()
    args.func(args)
