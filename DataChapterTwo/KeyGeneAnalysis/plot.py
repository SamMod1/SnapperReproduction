import pandas as pd
from pylab import plot, show, savefig, xlim, figure, \
                ylim, legend, boxplot, setp, axes, axvline
from datetime import datetime
from matplotlib.patches import Patch
import numpy as np
from src.stats import stats_test
from scipy import stats as scistat
from brokenaxes import brokenaxes
from matplotlib.pyplot import text
import string
import random
import statsmodels.stats.multicomp as multi
import scikit_posthocs as sp
from itertools import combinations


ALL_PS = [1.1928708820153783e-05, 1.19e-05, 2.46e-08, 5.81e-05, 4.19e-06,
          0.000428, 9.39e-06, 0.000158, 0.000136, 0.00225,
          0.308, 2.23e-05, 0.021359868498240783, 0.9059341426274077,
          0.4748140958549969, 9.485239525893153e-05, 6.5e-07, 0.000452,
          0.00271, 2.75e-07, 0.572, 7.45e-07, 0.703,
          1.99e-10, 7.13e-10, 4.83e-09, 0.00829628588516882, 2.5669719456581822e-06]

P_START = 0
EXPRESSION_FILE = "C:/Users/cfnsjm/Local Doccuments/PhD/Programming/R/DESeq_final/NormalisedData/normalised_counts_gonad_length_scaled.tsv"
GROUPS = ['F1', 'F3', 'F4', 'F5', 'F2']
#GROUPS = ['M3', 'M4', 'M2']
COLOURS = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']
GENES = ['star', 'cyp11a1.2', 'cyp17a1', 'hsd3b1', 'hsd17b1']
TRANSFORMATIONS = [None, np.log, np.log, np.log, np.log]
#TRANSFORMATIONS = [np.log, np.log, np.log, np.log, np.log]
#GENES = ['cyp19a1a', 'foxl2a', 'cyp11b1', 'hsd11b2']
#TRANSFORMATIONS = [np.log, np.log, None, None]
#TRANSFORMATIONS = [None, None, np.log, lambda x: x]
#GENES = ['fshr', 'lhcgr_ChrAur1', 'ar', 'armt1_fusion_2', 'esr2a']
#TRANSFORMATIONS = [np.log, np.log, np.log, np.log, lambda x: x]
#TRANSFORMATIONS = [np.log, np.log, np.log, None, None]
GENE_DISPLAY_NAMES = None
#GENE_DISPLAY_NAMES = ['fshr', 'lhcgr', 'ar', 'esr1', 'esr2a']
TEXT_SIZE = 10
Y_LIMS = (-0.2, 14.7)
LINE_WIDTH = .8
ALPHA = 0.05
ANNOTATE_SIGNIF = True
ALPHA_SIGNIF_LEVELS = {'0.05': '*', '0.01': '**', '0.001': '***'}
TRANSFORM = lambda x: np.log(x + 1)# lambda x: x
SAMPLES_TO_EXCULDE = ['X90_M3']


def multi_comparisons_letter_df_generator(comparisons_df, letter_ordering_series = None, 
                                          primary_optimisation_parameter = "Number of different letters", 
                                          monte_carlo_cycles = 5, letter_separator = '', ): 
    """
    Function takes a df listing pairwise comparisons with a cols labelled 'group1' and 'group2' for the two groups being compared 
    and another column labelled 'reject' with boolean values corresponding to whether the null hypothesis should be rejected 
    i.e. True: Both treatments are significantly different
    
    letter_ordering_series (default = None): In theory, which letters are assigned to each non-significance grouping is
    arbitrary and therefor the order can be changed. Offering letter_ordering_series a series with the same index as the output
    will make sure that the order that letters are assigned will follow letter_ordering_series from max to min. For boxplots,
    letter_ordering_series with median values is a good choice.
    
    monte_carlo_cycles (default = 5): Function will always return correct letter representation however it may be suboptimal. 
    Within each monte carlo cycle, a random letter is deleted until the representation breaks. The result with the optimum
    number layout of letters after n monte_carlo_cycles is returned. 
    
    The optimum letter layout is set by primary_optimisation_parameter (default = "Number of different letters"):
        'Number of different letter' optimises for fewest different letters
        "Min letters per row" optimises for the fewest letters assigned per treatment
        "Letter total" optimises for the fewest total letters of the treatments combined
        
    letter_separator (default = ''): Separator for each letter in string assigned to each treatment
    
    
    Letter representation is determined by the method described by Piepho 2004: An Algorithm for a Letter-Based Representation
    of All-Pairwise Comparisons
    """
    #'Insert' stage
    #make df with all unique groups as index
    letters_df = pd.concat((comparisons_df['group1'], comparisons_df['group2'])).drop_duplicates().to_frame().set_index(0)#comparisons_df['group1'].append(comparisons_df['group2']).drop_duplicates().to_frame().set_index(0)
    
    letters_df[letters_df.shape[1]] = 1
    for pos_result in comparisons_df.loc[comparisons_df['reject']==True].index:
        group1 = comparisons_df.loc[pos_result, 'group1']
        group2 = comparisons_df.loc[pos_result, 'group2']
        for letter_col in letters_df:
            group1_val = letters_df.loc[group1,letter_col]
            group2_val = letters_df.loc[group2,letter_col]
            if group1_val == 1 and group2_val == 1:
                #duplicate column
                new_col = letters_df.shape[1]
                letters_df[new_col] = letters_df[letter_col]
                #del val at group1 first col and at group2 new col
                letters_df.loc[group1,letter_col] = 0
                letters_df.loc[group2,new_col] = 0
    #'Absorb' stage          
    for col in letters_df:
       other_cols_list = list(letters_df)
       other_cols_list.remove(col)
       col_total = letters_df[col].sum()
       for other_col in other_cols_list:
           matched_total = 0
           for row in letters_df.index:
               if letters_df.loc[row, col] == 1 and letters_df.loc[row, other_col]: matched_total +=1
           if col_total == matched_total:
               letters_df.drop(col, axis = 1, inplace = True)  
               break
        
    def check_letters_against_tests(test_df, letters_df):
        if letters_df.sum(axis = 1).min() == 0: return False
        for result_row in test_df.index:
            group1 = test_df.loc[result_row, 'group1']
            group2 = test_df.loc[result_row, 'group2']
            reject = bool(test_df.loc[result_row, 'reject'])
            count_of_true_trues = 0
            count_of_true_falses = 0
            for letter_col in letters_df:
                group1_val = letters_df.loc[group1,letter_col]
                group2_val = letters_df.loc[group2,letter_col]
                if reject:
                    if group1_val != group2_val: count_of_true_trues += 1
                    if group1_val == 1 and group2_val == 1: 
                        return False
                if reject == False:
                    if group1_val == 1 and group2_val == 1: count_of_true_falses += 1
            if reject and count_of_true_trues == 0: 
                return False
            if reject == False and count_of_true_falses == 0: 
                return False
        return True

    #'Sweep stage' with monte carlo optimisation
    for i in range(monte_carlo_cycles):
        num_of_letters = letters_df.sum().sum()
        num_list = list(np.arange(start = 1, stop = 1+ num_of_letters))
        letters_df_monte_order = letters_df.copy()
        for row in letters_df_monte_order.index:
            for col in letters_df_monte_order:
                if letters_df_monte_order.loc[row,col] == 0: continue
                random_num = random.sample(num_list, 1)[0]
                letters_df_monte_order.loc[row,col] = random_num
                num_list.remove(random_num)
        
        current_letters_df = letters_df.copy()
        for pos in range(num_of_letters + 1):     
            mask = letters_df_monte_order.isin([pos])
            zero_df = letters_df.copy().loc[:] = 0
            letters_df_copy = current_letters_df.copy()
            letters_df_copy.mask(mask, other = zero_df, inplace = True)
            if check_letters_against_tests(comparisons_df,letters_df_copy):
                current_letters_df = letters_df_copy
        
        for col in letters_df:
            if current_letters_df[col].sum() == 0: current_letters_df.drop(col, axis = 1, inplace = True)
            
        # determine fitness parameters for optimisation
        current_fitness_parameter_vals = {"Min letters per row":current_letters_df.sum(axis = 1).max(),
                                          "Number of different letters": current_letters_df.shape[1],
                                          "Letter total": current_letters_df.sum().sum()}
        if i == 0: 
            best_fitness_parameter_vals = current_fitness_parameter_vals
            best_letters_df = current_letters_df
            continue
        
        if current_fitness_parameter_vals[primary_optimisation_parameter] > best_fitness_parameter_vals[primary_optimisation_parameter]:
            continue
        if current_fitness_parameter_vals[primary_optimisation_parameter] < best_fitness_parameter_vals[primary_optimisation_parameter]:
            best_letters_df = current_letters_df.copy()
            best_fitness_parameter_vals = current_fitness_parameter_vals
            
        if sum(current_fitness_parameter_vals.values()) < sum(best_fitness_parameter_vals.values()):
            best_letters_df = current_letters_df.copy()
            best_fitness_parameter_vals = current_fitness_parameter_vals
    
    #order cols
    if isinstance(letter_ordering_series, pd.Series):
        scoring_df = pd.DataFrame(index = best_letters_df.index)
        for row in best_letters_df.index:
            for col in best_letters_df:
                scoring_df.loc[row, col] = best_letters_df.loc[row, col] * letter_ordering_series[row]
        scoring_df = scoring_df.replace(0, np.NaN)
        scoring_means = scoring_df.mean(axis = 0).sort_values(ascending = False)
        best_letters_df = best_letters_df[scoring_means.index]
    # letter the cols     
    for col_name, col_num in zip(best_letters_df, range(len(best_letters_df.columns))):
        letter = string.ascii_lowercase[col_num]
        best_letters_df.loc[best_letters_df[col_name] == 1, col_name] = letter
    # make df with strings ready for presentation
    best_string_df = pd.DataFrame(index = best_letters_df.index)
    best_string_df.loc[:,'string'] = ""
    for row in best_letters_df.index:
        for col in best_letters_df:
            if best_letters_df.loc[row, col] != 0:
                letter_string = best_string_df.loc[row, 'string']
                letter = best_letters_df.loc[row, col]
                if letter_string == "": best_string_df.loc[row, 'string'] = letter
                else: best_string_df.loc[row, 'string'] = letter_separator.join((letter_string, letter))
                
    return best_string_df

def stack_correlation_table(df):
    """
    Converts a dataframe correlation table to a stacked comparisons table
    """
    df = df.stack().to_frame()
    for row in df.index:
        if row[0] == row[1]: 
            df = df.drop(row)
            continue
        sorted_row = list(row)
        sorted_row.sort()
        df.loc[row,'A'], df.loc[row,'B'] = sorted_row[0], sorted_row[1]
    df = df.set_index(['A', 'B'], drop = True)
    df = df.loc[~df.index.duplicated(keep='first')]
    return df

def scikit_results_munger(results, alpha):
    results = stack_correlation_table(results)
    results.rename({0:'p'}, axis = 1, inplace = True)
    results.loc[results['p'] <= alpha, 'reject'] = True
    results.loc[results['p'] > alpha, 'reject'] = False
    for row in results.index:
        results.loc[row, 'group1'] = row[0]
        results.loc[row, 'group2'] = row[1]
    return results

def post_hoc_df(df, Y_col, X_col, posthoc = "tukey", alpha = 0.05):
    """
    Returns a df with pairwise comparisons with reject column calculated according to alpha
    
    TODO: Add more posthoc tests to this function
    """
    if posthoc == "Statsmodels_tukey":
        comp = multi.MultiComparison(df[Y_col], df['comb'])
        results = comp.tukeyhsd(alpha=alpha)
        results = pd.DataFrame(data=results._results_table.data[1:], columns=results._results_table.data[0])
    if posthoc == "dunn": results = scikit_results_munger(sp.posthoc_dunn(df, val_col= Y_col, group_col= X_col, p_adjust = 'holm'), alpha)
    if posthoc == "tukey": results = scikit_results_munger(sp.posthoc_tukey(df, val_col= Y_col, group_col= X_col), alpha)
    return results


def build_signif_symbol(alpha_signif_symbols, p):
    levels = sorted(list(alpha_signif_symbols.keys()), key = lambda x: float(x))
    for signif_level in levels:
        if p <= float(signif_level):
            return alpha_signif_symbols[signif_level]
    
    raise ValueError("Value is not statistically significant. The alpha level or alpha significance levels symbol mappings may be wrong")
    
    
def check_overlap(positions, y, feature_height):
    positions = list(positions)
    for position in positions:
        if position-feature_height <= y-feature_height <= position+feature_height:
            new_y = y + (feature_height * 2.3)
            y = check_overlap(positions, new_y, feature_height)
  
    return y


corrected_ps = list(scistat.false_discovery_control(ALL_PS))

genes = [x.replace('gene-', '') for x in GENES]
expression_file = EXPRESSION_FILE
df = pd.read_csv(expression_file, sep='\t')

# df['gene_name'] = df['gene_id'].str.replace('gene-', '')#df.index.str.replace('gene-', '')
# df = df.drop(('gene_id'), axis = 1)
# df.index = range(len(df))
# #df = df.drop('gene_id', axis = 1)

df = df[df.columns[~df.columns.isin(SAMPLES_TO_EXCULDE)]]
groups = {group: {'samples': []} for group in GROUPS}
samples = np.array(df.columns[2:], dtype = str)
for group in GROUPS:
    group_samples = samples[[True if group in x else False for x in samples]]
    groups[group]['samples'].append(group_samples)


fig = figure()
ax = axes()
ax.set_ylim(Y_LIMS[0], Y_LIMS[1])
all_positions = []


all_ps = []

for gene_idx, gene in enumerate(genes):
    expression = df.loc[df['gene_name'] == gene].squeeze()[1:]
    for group in groups:
        group_expression = expression.loc[np.isin(expression.index, groups[group]['samples'])]
        groups[group][gene] = np.array(group_expression, dtype = 'double')
    
    positions = np.arange(1, len(GROUPS) + 1) + gene_idx * (len(GROUPS) + 1)
    all_positions.append(positions)
    gene_data = [groups[group][gene] for group in GROUPS]
    if np.any([np.any(x) for x in gene_data]):
        if TRANSFORMATIONS:
            stats_results, signif_groups = stats_test(GROUPS, gene_data, ALPHA, transformation = TRANSFORMATIONS[gene_idx])
        else:
            stats_results, signif_groups = stats_test(GROUPS, gene_data, ALPHA)
        all_ps.append(stats_results)
    else:
        all_ps.append(1)
    
    if corrected_ps[P_START + gene_idx] < ALPHA and signif_groups and signif_groups[0] and ANNOTATE_SIGNIF and len(GROUPS) > 2:
        reject = [x[:2] for x in signif_groups]
        signif_df = pd.DataFrame(signif_groups, columns = ('group1', 'group2', 'p'))
        group_combinations = np.array(list(combinations(GROUPS, 2)))
        signif_df = pd.DataFrame(group_combinations, columns = ('group1', 'group2'))
        signif_df['reject'] = pd.Series(list(combinations(GROUPS, 2))).isin(reject)
        median_df = pd.Series([x.mean() for x in gene_data])
        median_df.index = GROUPS
        
        group_annotations = list(multi_comparisons_letter_df_generator(signif_df, letter_ordering_series = median_df)['string'])
    
    gene_data = [TRANSFORM(x) for x in gene_data]
    plt = boxplot(
        gene_data, 
        positions = positions, 
        widths = .6, 
        patch_artist=True, 
        flierprops = {'markersize': 3, 'marker': '.', 'markerfacecolor': COLOURS[0]}
    )
    
    if corrected_ps[P_START + gene_idx] < ALPHA and signif_groups and signif_groups[0] and ANNOTATE_SIGNIF and len(GROUPS) > 2:
        for group_idx in range(len(GROUPS)):
            text(positions[group_idx] - abs(1-1.2/len(group_annotations[group_idx])), 
                 max(gene_data[group_idx]) + .5, group_annotations[group_idx], fontdict = {'size': TEXT_SIZE * .8})
    
    if not gene_idx == len(genes) - 1:
        axvline(positions[-1] + 1, linestyle = 'dashed', color = 'grey', lw=LINE_WIDTH)
    n_groups = len(GROUPS)
    
    for group in range(n_groups):
        colour = COLOURS[group]
        setp(plt['boxes'][group], color=colour)
        plt['boxes'][group].set_facecolor(colour)
        setp(plt['caps'][group * 2], color=colour)
        setp(plt['caps'][group * 2 + 1], color=colour)
        setp(plt['whiskers'][group * 2], color=colour)
        setp(plt['whiskers'][group * 2 + 1], color=colour)
        setp(plt['fliers'][group], markeredgecolor=colour, markerfacecolor=colour)
        setp(plt['medians'][group], color='black')
    
    if len(GROUPS) == 2:
        if ANNOTATE_SIGNIF and corrected_ps[P_START + gene_idx] < ALPHA:
            max_val = max([max(x) for x in gene_data])
            scale = max_val / 10
            line_positions = set([])
            feature_height = (LINE_WIDTH / 3) * scale
            for (group1, group2, p) in signif_groups:
                pos1 = GROUPS.index(group1)
                pos2 = GROUPS.index(group2)
                x1, x2 = positions[pos1], positions[pos2]
                y = max(max(gene_data[pos1]), max(gene_data[pos2])) + .3 * scale
                colour = 'black'
                
                y = check_overlap(line_positions, y, feature_height)
               
                line_positions.add(y)
     
                ax.plot([x1, x1, x2, x2], [y, y+feature_height, y+feature_height, y], lw=LINE_WIDTH, c=colour)
                ax.text(
                    (x1+x2)*.5, y+feature_height - feature_height / 4, build_signif_symbol(ALPHA_SIGNIF_LEVELS, corrected_ps[P_START + gene_idx]), ha='center', va='bottom', 
                    color=colour, fontsize=TEXT_SIZE * .7
                )
                
ax.set_xticks([np.median(x) for x in all_positions])
if GENE_DISPLAY_NAMES:
    ax.set_xticklabels(GENE_DISPLAY_NAMES, fontsize = TEXT_SIZE)
else:
    ax.set_xticklabels(GENES, fontsize = TEXT_SIZE)
ax.set_ylabel('Log Normalised Gene Expression')
print(all_ps)
print(list(corrected_ps[P_START:P_START+len(GENES)]))

            
        
# ====================== GENE GROUPS ============================
#GENES = ['star', 'cyp11a1.2', 'cyp17a1', 'hsd3b1', 'hsd17b1']
#GENES = ['cyp19a1a', 'foxl2a', 'cyp11b1', 'hsd11b2']
        