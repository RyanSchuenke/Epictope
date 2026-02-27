import pandas as pd
import matplotlib.pyplot as plt
import os
import logging
logger = logging.getLogger(__name__)

def plot_scores(scores_file: os.PathLike, output_file:os.PathLike = None) -> None:
    """
    Plots the moving average of the min score determined by epictope
    
    :param scores_file: scores csv file produced by the epictope pipeline
    :type scores_file: os.PathLike
    :param output_file: output destination of the figure
    :type output_file: os.PathLike
    """
    if not os.path.exists(scores_file):
        logger.error("Score file '"+scores_file+"' not found")
        raise Exception("Score file '"+scores_file+"' not found")
    if not output_file:
        output_file = os.path.splitext(scores_file)[0]+".png"
    
    plot_df = pd.read_csv(scores_file).set_index(["position", "aa"])
    fig, ax = plt.subplots(1,1)
    fig.set_size_inches(12,6)
    ax.set_xlim((1, plot_df.index.get_level_values("position").max()))
    ax.set_ylim((-0.05, 1.05))
    ax.minorticks_on()
    ax.set_xlabel("Amino Acid Position")
    ax.set_ylabel("Score")
    
    plot_df = plot_df.dropna()
    window = 7
    half_window = window//2
    
    # get the moving average with a window of size 7 (3, 1, 3) around each position
    for i, row in plot_df[:half_window].iterrows():
        index = i[0]
        plot_df.loc[index,"min_moving_avg"] = plot_df.loc[:index+half_window-1]["min"].mean()
    for i, row in plot_df[half_window:len(plot_df)-half_window].iterrows():
        index = i[0]
        plot_df.loc[index,"min_moving_avg"] = plot_df.loc[index-half_window:index+half_window]["min"].mean()
    for i, row in plot_df[len(plot_df)-half_window:].iterrows():
        index = i[0]
        plot_df.loc[index,"min_moving_avg"] = plot_df.loc[index-half_window:]["min"].mean()
    
    ax.plot(plot_df.index.get_level_values("position"), plot_df["min_moving_avg"])
    fig.savefig(output_file)

def plot_score_components(scores_file: os.PathLike, output_file:os.PathLike = None, 
                          score_components:list[str] = ["rsa", "ss_score", "inv_anchor2", "normalized_entropy"]) -> None:
    """
    Plots the moving average of the four score components on the same plot, as determined by epictope
    
    :param scores_file: scores csv file produced by the epictope pipeline
    :type scores_file: os.PathLike
    :param output_file: output destination of the figure
    :type output_file: os.PathLike
    :param score_components: list of score component columns to plot
    :type score_components: list[str]
    """
    if not os.path.exists(scores_file):
        logger.error("Score file '"+scores_file+"' not found")
        raise Exception("Score file '"+scores_file+"' not found")
    if not output_file:
        output_file = os.path.splitext(scores_file)[0]+"_components.png"
    
    plot_df = pd.read_csv(scores_file).set_index(["position", "aa"])
    fig, ax = plt.subplots(1,1)
    fig.set_size_inches(12,6)
    ax.set_xlim((1, plot_df.index.get_level_values("position").max()))
    ax.set_ylim((-0.05, 1.05))
    ax.minorticks_on()
    ax.set_xlabel("Amino Acid Position")
    ax.set_ylabel("Score")
    
    plot_df = plot_df.dropna()
    window = 7
    half_window = window//2
    
    # get the moving average with a window of size 7 (3, 1, 3) around each position
    for col in score_components:
        for i, row in plot_df[:half_window].iterrows():
            index = i[0]
            plot_df.loc[index,col+"_moving_avg"] = plot_df.loc[:index+half_window-1][col].mean()
        for i, row in plot_df[half_window:len(plot_df)-half_window].iterrows():
            index = i[0]
            plot_df.loc[index,col+"_moving_avg"] = plot_df.loc[index-half_window:index+half_window][col].mean()
        for i, row in plot_df[len(plot_df)-half_window:].iterrows():
            index = i[0]
            plot_df.loc[index,col+"_moving_avg"] = plot_df.loc[index-half_window:][col].mean()
        
        ax.plot(plot_df.index.get_level_values("position"), plot_df[col+"_moving_avg"], label=col)
    fig.legend()
    fig.savefig(output_file)
