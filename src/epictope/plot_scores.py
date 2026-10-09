import pandas as pd
import matplotlib.pyplot as plt
import os
import logging
logger = logging.getLogger(__name__)

def plot_scores(scores_file: os.PathLike, output_file:os.PathLike = None, 
                score_components:list[str] = ["min"], 
                window_size: int = 7, min_non_zero_count: int = 1) -> None:
    """
    Plots the moving average of the four score components on the same plot, as determined by EpicTope
    
    :param scores_file: scores csv file produced by the EpicTope pipeline
    :type scores_file: os.PathLike
    :param output_file: output destination of the figure
    :type output_file: os.PathLike
    :param score_components: list of score component columns to plot
    :type score_components: list[str]
    :param window_size: size of the sliding window used to average the EpicTope min score
    :type window_size: int
    :param min_non_zero_count: the number of non-zero min scores needed for the average score to be calculated, otherwise set to zero
    :type min_non_zero_count: int
    """
    if not os.path.exists(scores_file):
        logger.error("Score file '"+scores_file+"' not found")
        raise Exception("Score file '"+scores_file+"' not found")
    if not output_file:
        output_file = os.path.splitext(scores_file)[0]+".png"
        out_csv = os.path.join(os.path.dirname(scores_file),os.path.splitext(os.path.basename(scores_file))[0]+"_plot_data.csv")
    else:
        out_csv = os.path.splitext(output_file)[0]+"_plot_data.csv"
    
    plot_df = pd.read_csv(scores_file).set_index(["position", "aa"])
    fig, ax = plt.subplots(1,1)
    fig.set_size_inches(12,6)
    ax.set_xlim((1, plot_df.index.get_level_values("position").max()))
    ax.set_ylim((-0.05, 1.05))
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.minorticks_on()
    ax.set_xlabel("Amino Acid Position", fontsize=18)
    ax.set_ylabel("Score", fontsize=18)
    
    plot_df = plot_df.dropna()
    half_window = window_size//2
    out_cols = []
    # get the moving average with a window of size 7 (3, 1, 3) around each position
    for col in score_components:
        for i, row in plot_df[:half_window].iterrows():
            index = i[0]
            plot_df.loc[index,col+"_moving_avg"] = plot_df.loc[:index+half_window-1][col].mean() if sum(1 for i in plot_df.loc[:index+half_window-1][col] if i) >= min_non_zero_count else 0
        for i, row in plot_df[half_window:len(plot_df)-half_window].iterrows():
            index = i[0]
            plot_df.loc[index,col+"_moving_avg"] = plot_df.loc[index-half_window:index+half_window][col].mean() if sum(1 for i in plot_df.loc[index-half_window:index+half_window][col] if i) >= min_non_zero_count else 0
        for i, row in plot_df[len(plot_df)-half_window:].iterrows():
            index = i[0]
            plot_df.loc[index,col+"_moving_avg"] = plot_df.loc[index-half_window:][col].mean() if sum(1 for i in plot_df.loc[index-half_window:][col] if i) >= min_non_zero_count else 0
        
        ax.plot(plot_df.index.get_level_values("position"), plot_df[col+"_moving_avg"], label=col)
        out_cols.append(col)
        out_cols.append(col+"_moving_avg")
    ax.legend(bbox_to_anchor=(0., .97, 1., .102), loc='upper center', ncols=len(score_components), mode="expand", borderaxespad=0., fontsize=12)
    
    plot_df[out_cols].to_csv(out_csv)
    fig.savefig(output_file, bbox_inches='tight')
