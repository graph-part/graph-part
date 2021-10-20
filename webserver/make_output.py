'''
Script to produce markdown output for nice
display in the GUI.
'''
import json
import pandas as pd
from tabulate import tabulate
import matplotlib.pyplot as plt

results = json.load(open("output/graphpart_report.json","r"))


def make_progress_plot(df, out_name):
    '''
    Plot removal progess as figure.
    Connectivity vs. iteration
    Number vs. iteration
    '''
    plt.figure()
    fig, ax1 = plt.subplots()

    color = 'tab:red'
    ax1.set_xlabel('Removal step')
    ax1.set_ylabel('Connectivity', color=color)
    ax1.plot(df.index, df['Connectivity'], color=color)
    ax1.tick_params(axis='y', labelcolor=color)

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    color = 'tab:blue'
    ax2.set_ylabel('Entitities', color=color)  # we already handled the x-label with ax1
    ax2.plot(df.index, df['#Entities'], color=color)
    ax2.tick_params(axis='y', labelcolor=color)

    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.savefig(out_name)


def make_partition_plot(df, out_name):
    df = df.set_index('label', drop=True).drop(['count'], axis=1)
    df.T.plot.bar(stacked=True)
    plt.xlabel('Partitions')
    plt.ylabel('Sequences')
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
make_partition_plot(df, 'partition_plot.png')
output_md += f'\n\n ![plot]({"partition_plot.png"})'
output_md += '\n\n'

## make detailed table of removal processes
if 'removal_step_1' in results:

    output_md += '### Homology removal report\n\n'
    df = pd.DataFrame.from_dict(results['removal_step_1']).T

    make_progress_plot(df, 'removal_step_1_plot.png')
    output_md += '#### Graph connectivity vs. number of sequences\n\n'
    output_md += f'\n\n ![plot]({"removal_step_1_plot.png"})'
    output_md += '\n\n'
    output_md += '#### Full removal metrics\n\n'
    output_md += tabulate(df, tablefmt='github', headers='keys')
    output_md += '\n\n'

if 'removal_step_2' in results:
    output_md += '### Homology removal report - second pass'
    df = pd.DataFrame.from_dict(results['removal_step_2']).T

    make_progress_plot(df, 'removal_step_2_plot.png')
    output_md += f'\n\n ![plot]({"removal_step_2_plot.png"})'
    output_md += '\n\n'
    output_md += tabulate(df, tablefmt='github', headers='keys')
    output_md += '\n\n'







open("output.md", "w").write(output_md)