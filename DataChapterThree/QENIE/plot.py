import pandas as pd
import numpy as np
import stats
import matplotlib.pyplot as plt
import seaborn as sns


SEED = 1
SUMMARY_FUNC = np.median
ALPHA = 0.05
FIGSIZE = [16, 6]
FONTSIZE = 24
YLIMS = (0.78, 1.22)
SEX_COLOURS = {'M': 'blue', 'F': 'red', 'J': 'grey', 'T': 'orange'}
ORDER = (
    'brain_to_hypo', 'brain_to_pit', 'brain_to_gonad', 'hypo_to_brain', 'hypo_to_pit', 'hypo_to_gonad',
    'pit_to_brain', 'pit_to_hypo', 'pit_to_gonad', 'gonad_to_brain', 'gonad_to_hypo', 'gonad_to_pit'
)
READABLE_LABELS = (
    'brain to hypothalamus', 'brain to pituitary', 'brain to gonad',
    'hypothalamus to brain', 'hypothalamus to pituitary', 'hypothalamus to gonad',
    'pituitary to brain', 'pituitary to hypothalamus', 'pituitary to gonad',
    'gonad to brain', 'gonad to hypothalamus', 'gonad to pituitary'
)
READABLE_LABELS = (
    'B -> H', 'B -> P', 'B -> G', 'H -> B', 'H -> P', 'H -> G',
    'P -> B', 'P -> H', 'P -> G', 'G -> B', 'G -> H', 'G -> P'
)
ALPHA_SIGNIF_LEVELS = {'0.05': '*', '0.01': '**', '0.001': '***'}


np.random.seed(SEED)

avg_data = 'DataChapterThree/QENIE/Results/All Sex/AverageSsec.csv'
combined_all = 'DataChapterThree/QENIE/Results/All Sex/CombinedSsec.csv'
male_all_data = 'DataChapterThree/QENIE/Results/Males/CombinedSsec.csv'
male_avg_data = 'DataChapterThree/QENIE/Results/Males/AverageSsec.csv'
female_all_data = 'DataChapterThree/QENIE/Results/Females/CombinedSsec.csv'
female_avg_data = 'DataChapterThree/QENIE/Results/Females/AverageSsec.csv'

avg = pd.read_csv(avg_data)
comb = pd.read_csv(combined_all)

all_m = pd.read_csv(male_all_data)
avg_m = pd.read_csv(male_avg_data)

all_f = pd.read_csv(female_all_data)
avg_f = pd.read_csv(female_avg_data)

combined_sexes = pd.concat([all_f, all_m])


def build_signif_symbol(alpha_signif_symbols, p):
    levels = sorted(list(alpha_signif_symbols.keys()), key = lambda x: float(x))
    for signif_level in levels:
        if p <= float(signif_level):
            return alpha_signif_symbols[signif_level]
    
    raise ValueError("Value is not statistically significant. The alpha level or alpha significance levels symbol mappings may be wrong")


def bootstrap_median_ci(data, n_bootstraps=5000, ci=95):
    """Calculate bootstrapped confidence intervals for the median."""
    medians = []
    for _ in range(n_bootstraps):
        # Resample with replacement
        resampled = np.random.choice(data, size=len(data), replace=True)
        medians.append(np.median(resampled))
    
    lower = np.percentile(medians, (100 - ci) / 2)
    upper = np.percentile(medians, ci + (100 - ci) / 2)
    return lower, upper


def nsec_plot(
        nsec_df, 
        sex_and_colours, 
        summary_func = np.mean,
        order = None,
        ylims = None,
        measurement_to_plot = 'value',
        alpha = 0.05,
        alpha_signif_levels = {'0.05': '*', '0.01': '**', '0.001': '***'},
        x_label = 'Tissue interaction', 
        y_label = r"$nS_{sec}$",
        plot_title = '',
        fontsize = 18,
    ):
    if not order:
        order = sorted(nsec_df['tissue'].unique())
    
    nsec_df['Sex'] = nsec_df['tissue'].str[-1]
    nsec_df['tissue'] = nsec_df['tissue'].str[:-1]
    sexes = nsec_df['Sex'].unique()
    
    # Preparing the data

    nsec_df = nsec_df.groupby(['tissue', 'Sex'])[measurement_to_plot].apply(list).reset_index()

    # print('=================================================')
    # print(results)
    # print('=================================================')
    
    # Unique months
    unique_tissues = sorted(nsec_df['tissue'].unique(), key = lambda x: order.index(x))

    # Box plot positions and widths
    box_width = 0.9  # Width of each box
    month_gap = 0.3  # Gap between groups of months
    position = 0
    starting_position = 0.3 # Starting position

    all_positions = []  # Collect all positions for setting x-ticks later

    # Plot setup
    fig, ax1 = plt.subplots(figsize=(14, 7))
    
    if ylims:
        plt.ylim(ylims)
    
    tissue_means = {}
    for tissue in unique_tissues:
        tissue_data = nsec_df[nsec_df['tissue'] == tissue]
        sex_means = []
        sex_positions = []
        sex_errors = []
        sex_data = []
        for sex_index, sex in enumerate(sexes):
            data = tissue_data[tissue_data['Sex'] == sex][measurement_to_plot].tolist()
            if data:  # Check if there's data to plot
                data[0] = [float(x) for x in data[0]]
                sex_data.append(data[0])
                mean = summary_func(data[0])
                sex_means.append(mean)
                tissue_means[sex] = tissue_means.get(sex, []) + [mean]
                quarts = bootstrap_median_ci(data[0], ci = (1-alpha) * 100)
                this_position = position + starting_position + sex_index * box_width
                ax1.plot(
                    this_position,
                    mean,
                    color=SEX_COLOURS[sex], 
                    marker='o'
                )
                sex_positions.append(this_position)
                ax1.errorbar(position + starting_position + sex_index * box_width, mean, yerr=((abs(mean - quarts[0]),),(abs(mean - quarts[1]),)), ecolor=SEX_COLOURS[sex], capsize=7, capthick=3, linewidth = 3)
                sex_errors.append(quarts[1])
                #y_data = tissue_data[tissue_data['Sex'] == sex][measurement_to_plot].explode().astype(float)
                #x_data = np.full(y_data.shape, position + sex_index * box_width) + starting_position
                #ax1.scatter(x_data, y_data, alpha=1, color=SEX_COLOURS[sex], edgecolor='none', s=20, zorder=3)
        ax1.plot(sex_positions, sex_means, color='black', linestyle='--', linewidth = 3)
        stat, p, d_of_f = stats.stats_test_two_groups(sex_data[0], sex_data[1], alpha)#results_dict.get(tissue + "M", {}).get(tissue + "F", 1)
        if p <= alpha:
            ax1.text(
                position + 0.7, max(sex_errors) + max(sex_errors) / 50, build_signif_symbol(alpha_signif_levels, p), ha='center', va='bottom', 
                color='black', fontsize=fontsize
            )
        
        position += box_width * len(sexes) + month_gap
        all_positions.append((position - (box_width * len(sexes) + month_gap) / 2) - month_gap)
    # Setting x-ticks to be in the middle of each group
    ax1.set_xticks(all_positions)
    ax1.set_xticklabels([x for x in READABLE_LABELS], rotation=45)
    # Final plot adjustments
    ax1.set_xlabel(x_label, size = fontsize)
    ax1.tick_params(axis='x', labelsize = fontsize * 0.8)
    ax1.set_ylabel(y_label, size = fontsize)
    ax1.tick_params(axis='y', labelsize = fontsize * 0.8)

    plt.tight_layout()
    plt.show()
    

combined_sexes['tissue'] = combined_sexes['tissue'].str.replace('mid', 'hypo')


def build_heatmap(sex, results_dict, unique_tissues):
    heatmap_data = []
    
    for tissue1 in unique_tissues:
        row_f = []
        row_m = []

        for tissue2 in unique_tissues:
            try:
                row_f.append(results_dict[tissue1 + sex][tissue2 + 'F'])
            except KeyError:
                row_f.append(1)
            try:
                row_m.append(results_dict[tissue1 + sex][tissue2 + 'M'])
            except KeyError:
                row_m.append(1)
                
        row = row_f + row_m
        heatmap_data.append(row)
        
    return heatmap_data


nsec_df = combined_sexes.copy(deep=True)
order = ORDER
measurement_to_plot = 'value'
alpha = ALPHA
readable_labels = READABLE_LABELS


if not order:
    order = sorted(nsec_df['tissue'].unique())

if not readable_labels:
    readable_labels = order

nsec_df['Sex'] = nsec_df['tissue'].str[-1]
nsec_df['tissue'] = nsec_df['tissue'].str[:-1]
sexes = nsec_df['Sex'].unique()

# Preparing the data

nsec_df = nsec_df.groupby(['tissue', 'Sex'])[measurement_to_plot].apply(list).reset_index()

results = stats.stats_test(nsec_df['tissue'] + nsec_df['Sex'], nsec_df[measurement_to_plot], alpha)

unique_tissues = sorted(nsec_df['tissue'].unique(), key = lambda x: order.index(x))

results_dict = {}
for result in results[1]:
    try:
        results_dict[result[0]][result[1]] = result[2]
    except KeyError:
        results_dict[result[0]] = {result[1]: result[2]}
    try:
        results_dict[result[1]][result[0]] = result[2]
    except KeyError:
        results_dict[result[1]] = {result[0]: result[2]}

heatmap_data = build_heatmap('F', results_dict, unique_tissues)
heatmap_data += build_heatmap('M', results_dict, unique_tissues)
heatmap_data = np.array(heatmap_data)
mask = np.triu(heatmap_data)
heatmap_data = heatmap_data <= alpha

labels = [x + ' (F)' for x in readable_labels] + [x + ' (M)' for x in readable_labels]
heatmap_data = pd.DataFrame(heatmap_data, columns = labels)
heatmap_data.index = labels

sns.heatmap(~heatmap_data, mask = mask)


#nsec_plot(combined_sexes.copy(deep=True), SEX_COLOURS, order = ORDER, fontsize = FONTSIZE)
nsec_plot(combined_sexes.copy(deep=True), SEX_COLOURS, SUMMARY_FUNC, ORDER, ylims = YLIMS, fontsize = FONTSIZE)
