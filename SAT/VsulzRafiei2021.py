import os
from myFunctions import load_data, reduce_word, merge_pdf_files
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from os.path import join as join
from scipy import stats as ss
# Switch under the working directory
os.chdir(r'/Users/bs3667/Dropbox (NYU Langone Health)/Bo Shen Working files/Speed-accuracy-tradeoff/Rafiei_Rahnev2021')
# Define I/O directories
datadir = r'/Users/bs3667/Dropbox (NYU Langone Health)/Bo Shen Working files/Speed-accuracy-tradeoff/Rafiei_Rahnev2021/data'
svdir = r'/Users/bs3667/Dropbox (NYU Langone Health)/Bo Shen Working files/Speed-accuracy-tradeoff/Rafiei_Rahnev2021/Analysis'
# loading data
mt = pd.read_csv(join(datadir, 'Transfromed.txt'), sep=',')
# See RT skewness
sublist = np.unique(mt.Subject)
dtype = [('subID', int), ('FisherPearson', float), ('TimePressure', int)]  # Adjust data types as needed
Skewness = np.empty(0, dtype=dtype)
subj = 0
pdf_files = []
sns.set_palette('Paired')
while subj < len(sublist):
    if subj % 12 == 0:
        fig, axs = plt.subplots(4, 3, figsize=(8.5, 11))
        lgd = True
    for r in range(4):
        for c in range(3):
            for i in range(5):
                RT = mt[(mt['Subject']==sublist[subj]) & (mt['condition'] == i+1)]['rt']
                data_to_append = (sublist[subj], ss.skew(RT, axis=None, bias=True, nan_policy='omit'), i+1)
                Skewness = np.append(Skewness, np.array(data_to_append, dtype=dtype))
                num_bins = 50
                bins = np.linspace(RT.min(), RT.max(), num_bins)
                sns.histplot(RT, bins=bins, kde=True, alpha=0.5, label=f'{i+1}', ax=axs[r, c])
            if lgd:
                axs[r,c].legend(fontsize='5')
            lgd=False
            axs[r,c].set_title(f'#{subj+1}')
            subj += 1
            if subj >= len(sublist):
                break
        else:
            continue
        break
    plt.tight_layout()
    filename = join(svdir, 'IndvRTSkewness_toSubj' + str(subj) + '.pdf')
    pdf_files.append(filename)
    plt.savefig(filename, format='pdf')
    plt.close(fig)
output_file = join(svdir, 'IndvRTSkewness_Cmbnd.pdf')
merge_pdf_files(pdf_files, output_file)
df = pd.DataFrame(Skewness)
df.to_csv(join(svdir, 'IndvRTSkewness.csv'), index=False)

# plotting
df = pd.read_csv(join(svdir, 'IndvRTSkewness.csv'))
wide_df = df.pivot(index='subID', columns='TimePressure', values='FisherPearson')
# Reset column names
wide_df.columns.name = None
wide_df.reset_index(inplace=True, drop=True)
wide_df = wide_df[wide_df.columns[::-1]]
# pair plot contrasting each condition pairs
def plot_unity(xdata, ydata, **kwargs):
    mn = -1
    mx = 25
    points = np.linspace(mn, mx, 100)
    plt.gca().plot(points, points, marker=None,
            linestyle='--', color='gray', label='Diagonal Line', zorder=1)
    plt.gca().set_xlim((mn, mx))
    plt.gca().set_ylim((mn, mx))
def hide_current_axis(*args, **kwds):
    plt.gca().set_visible(False)
sns.set_palette('Paired')
g = sns.pairplot(wide_df, corner=True, plot_kws=dict(alpha=0.8, s=20))
g.map_lower(plot_unity)
g.map_diag(hide_current_axis)
for i in range(len(g.axes)):
    g.axes[i, i].axis('off')
g.fig.set_size_inches(7,7)
plt.tight_layout()
filename = join(svdir, 'IndvRTSkewness.pdf')
plt.savefig(filename, format='pdf')

# Initialize a grid of plots with an Axes for each walk
grid = sns.FacetGrid(df, col="subID", hue="subID", palette="tab20c",
                     col_wrap=5, height=1.5)
# Draw a horizontal line to show the starting point
grid.refline(y=0, linestyle=":")
# Draw a line plot to show the trajectory of each random walk
grid.map(plt.plot, "TimePressure", "FisherPearson", marker="o")
# Adjust the tick positions and labels
grid.set(xticks=np.arange(1,6), xlabel='Condition', ylabel='Skewness',
         xlim=(0, 5), ylim=(-1, 25))
# Adjust the arrangement of the plots
grid.fig.tight_layout(w_pad=1)
filename = join(svdir, 'IndvRTSkewness_cond.pdf')
plt.savefig(filename, format='pdf')
plt.close(fig)

fig = plt.figure(figsize=(3, 3))
sns.set_palette('Paired')
plt.plot([min(df['FisherPearson']), max(df['FisherPearson'])],
         [min(df['FisherPearson']), max(df['FisherPearson'])],
         color='gray', linestyle='--', label='Diagonal Line', zorder=1)
plt.scatter(df[df['TimePressure'] == 'Low']['FisherPearson'],
            df[df['TimePressure'] == 'High']['FisherPearson'], alpha=.8, s=7)
plt.title('RT skewness (Fisher-Pearson)')
plt.xlabel('Low time pressure')
plt.ylabel('High time pressure')
plt.tight_layout()
filename = join(svdir, 'Individual_choice', 'IndvRTSkewness.pdf')
plt.savefig(filename, format='pdf')
plt.close(fig)
