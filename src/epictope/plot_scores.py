import pandas as pd
import matplotlib.pyplot as plt
import os
import logging
logger = logging.getLogger(__name__)
plt.style.use('tableau-colorblind10')

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
    
    # initialize transmembrane state
    states = {'i': "Inside membrane",
              'm': "Transmembrane",
              'o': "Outside membrane",
              'n': "Signal peptide",
              'h': "Signal peptide",
              'c': "Signal peptide",
              '/': "Signal peptide"}
    states_colors = {"Inside membrane":"#5F9ED1",
                     "Transmembrane":"#898989",
                     "Outside membrane":"#A2C8EC",
                     "Signal peptide":"#FFBC79"}
    membrane_state_start = 1
    membrane_state = states[plot_df.loc[1]["phobius"].values[0].lower()]
    membrane_handles = []
    
    # get the moving average with a window of size 7 (3, 1, 3) around each position
    for index in range(1,len(plot_df)+1):
        # get moving average of scores to graph
        if index <= half_window:
            for col in score_components:
                plot_df.loc[index,col+"_moving_avg"] = plot_df.loc[:index+half_window-1][col].mean() if sum(1 for i in plot_df.loc[:index+half_window-1][col] if i) >= min_non_zero_count else 0
        elif index <= len(plot_df)-half_window:
            for col in score_components:
                plot_df.loc[index,col+"_moving_avg"] = plot_df.loc[index-half_window:index+half_window][col].mean() if sum(1 for i in plot_df.loc[index-half_window:index+half_window][col] if i) >= min_non_zero_count else 0
        elif index <= len(plot_df):
            for col in score_components:
                plot_df.loc[index,col+"_moving_avg"] = plot_df.loc[index-half_window:][col].mean() if sum(1 for i in plot_df.loc[index-half_window:][col] if i) >= min_non_zero_count else 0
        # get transmembrane regions to highlight
        if states[plot_df.loc[index]["phobius"].values[0].lower()] != membrane_state:
            membrane_handles.append(ax.axvspan(membrane_state_start, index, facecolor=states_colors[membrane_state], alpha=0.4, label=membrane_state))
            membrane_state = states[plot_df.loc[index]["phobius"].values[0].lower()]
            membrane_state_start = index
    if membrane_state_start != 1:
        membrane_handles.append(ax.axvspan(membrane_state_start, index, facecolor=states_colors[membrane_state], alpha=0.4, label=membrane_state))
    # plot moving averages
    score_handles = []
    for col in score_components:
        line, = ax.plot(plot_df.index.get_level_values("position"), plot_df[col+"_moving_avg"], label=col)
        score_handles.append(line)
        out_cols.append(col)
        out_cols.append(col+"_moving_avg")
    
    if score_handles:
        ax.legend(score_handles, score_components,
                                 bbox_to_anchor=(0., 1.08 if membrane_handles else .97, 1., .102),
                                 loc='upper center', ncols=len(score_components),
                                 mode="expand", borderaxespad=0., fontsize=12)
    if membrane_handles:
        ax2 = ax.twinx()
        ax2.get_yaxis().set_ticks([])
        membrane_labels = [handle.get_label() for handle in membrane_handles]
        by_label = dict(zip(membrane_labels, membrane_handles))
        ax2.legend(by_label.values(), by_label.keys(),
                                    bbox_to_anchor=(0., .97, 1., .102), loc='upper center',
                                    ncols=len(by_label), mode="expand", borderaxespad=0., fontsize=12)
    # save output image and data
    plot_df[out_cols].to_csv(out_csv)
    fig.savefig(output_file, bbox_inches='tight')
