'''
Script to produce markdown output for nice
display in the GUI.
'''
import json
import pandas as pd
from tabulate import tabulate
import matplotlib.pyplot as plt
import datetime

results = json.load(open("output/graphpart_report.json","r"))


def make_progress_plot(df, threshold, out_name):
    '''
    Plot removal progess as figure.
    Connectivity vs. iteration
    Number vs. iteration
    '''
    fig, (fig1_ax1, fig2_ax1) = plt.subplots(1,2, figsize = (12,4.5))

    color = 'tab:red'
    fig1_ax1.set_xlabel('Removal step')
    fig1_ax1.set_ylabel('Connectivity', color=color)
    fig1_ax1.plot(df.index, df['Connectivity'], color=color)
    fig1_ax1.tick_params(axis='y', labelcolor=color)

    ax2 = fig1_ax1.twinx()  # instantiate a second axes that shares the same x-axis

    color = 'tab:blue'
    ax2.set_ylabel('Remaining sequences', color=color)  # we already handled the x-label with ax1
    ax2.plot(df.index, df['#Entities'], color=color)
    ax2.tick_params(axis='y', labelcolor=color)


    color = 'tab:red'
    fig2_ax1.set_xlabel('Removal step')
    fig2_ax1.set_ylabel('# Problematic sequences', color=color)
    fig2_ax1.plot(df.index, df['#Problematics'], color=color)
    fig2_ax1.tick_params(axis='y', labelcolor=color)

    ax2 = fig2_ax1.twinx()  # instantiate a second axes that shares the same x-axis

    color = 'tab:blue'
    ax2.set_ylabel('Minimal distance between partitions', color=color)  # we already handled the x-label with ax1
    ax2.plot(df.index, df['Min-threshold'], color=color)
    ax2.tick_params(axis='y', labelcolor=color)
    ax2.axhline(threshold, linestyle='--', label='threshold (transformed)')
    ax2.legend()
    
    plt.suptitle('Iterative sequence removal procedure')
    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.savefig(out_name)


def make_partition_plot(df, out_name):
    df = df.set_index('label', drop=True).drop(['count'], axis=1)
    df.T.plot.bar(stacked=True, figsize = (12,4.5))
    plt.xlabel('Partitions')
    plt.ylabel('Sequences')
    plt.title('Composition of partitions')
    plt.tight_layout()
    plt.savefig(out_name)


output_md = '' # collect all output elements in this string.

output_md += f'## Graph-Part - Results\n\n'

output_md += f'[Download partition assignments](graphpart_result.csv)\n\n'

df = pd.DataFrame.from_dict(results['labels_start']).T
output_md += '### Summary of input data\n\n'
output_md += f'- Received **{results["samples_pre_removal"]:,}** sequences with **{len(df)}** classes.\n\n'
output_md += f'- Found **{results["graph_edges_start"]:,}** pairwise identities higher than the threshold **{results["config"]["threshold"]}**.\n\n'
output_md += tabulate(df, tablefmt="github", headers=['Label', 'i', 'n', 'n_target'])
output_md += '\n\n'

output_md += f'### Assignment summary\n\n'
output_md += '\n\n'
output_md += f'Retained **{results["samples_after_removal"]:,}** sequences.\n\n'

df = pd.read_json(results['partitioning_after_removal'])
output_md += tabulate(df, tablefmt='github', headers='keys')
output_md += '\n\n'
make_partition_plot(df, 'output/partition_plot.png')
output_md += f'\n\n ![plot]({"partition_plot.png"})'
output_md += '\n\n'



delta = datetime.timedelta(seconds=results['time_script_complete'] - results['time_script_start'])
output_md += f'Total runtime: {str(delta)}\n\n'

delta = datetime.timedelta(seconds=results['time_edges_complete'] - results['time_script_start'])
output_md += f'Alignment runtime: {str(delta)}\n\n'






## make detailed table of removal processes
if 'removal_step_1' in results:

    output_md += '### Homology removal report\n\n'
    df = pd.DataFrame.from_dict(results['removal_step_1']).T

    make_progress_plot(df, results['threshold_transformed'], 'output/removal_step_1_plot.png')
    output_md += '#### Graph connectivity vs. number of sequences\n\n'
    output_md += f'\n\n ![plot]({"removal_step_1_plot.png"})'
    output_md += '\n\n'
    output_md += '#### Full removal metrics\n\n'
    output_md += tabulate(df, tablefmt='github', headers='keys')
    output_md += '\n\n'

if 'removal_step_2' in results:
    output_md += '### Homology removal report - second pass'
    df = pd.DataFrame.from_dict(results['removal_step_2']).T

    make_progress_plot(df, results['threshold_transformed'], 'output/removal_step_2_plot.png')
    output_md += f'\n\n ![plot]({"removal_step_2_plot.png"})'
    output_md += '\n\n'
    output_md += tabulate(df, tablefmt='github', headers='keys')
    output_md += '\n\n'







open("output.md", "w").write(output_md)