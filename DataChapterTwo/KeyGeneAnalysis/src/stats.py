import pandas as pd
import numpy as np
import re
import scikit_posthocs as sp
import scipy.stats as stats
import math
from itertools import combinations
from statsmodels.stats.multicomp import pairwise_tukeyhsd


TISSUE = "gonad"
SAMPLE_NAME_REGEX = re.compile(
    "[0-9]+_[A-z]+\.([0-9]+|[A-z]+)_[A-z]+\.([A-z]|[0-9])"
)  # Regex used to identify sample names
REGEX_SPLIT = re.compile(
    "^X[0-9]+_"
)  # Regex used to remove the sample number identifier
INCLUDE = [
    "SEX.FEMALE_STAGE.1",
    "SEX.FEMALE_STAGE.3",
    "SEX.FEMALE_STAGE.5",
    "SEX.FEMALE_STAGE.6",
]
TRANSFORMATIONS = [
    math.log,
    math.log10,
    math.sqrt,
    lambda x: 1 / x,
    lambda x: math.asin(math.sqrt(x)),
]
LETTERS = [
    "a",
    "b",
    "c",
    "d",
    "e",
    "f",
    "g",
    "h",
    "i",
    "j",
    "k",
    "l",
    "m",
    "n",
    "o",
    "p",
    "q",
    "r",
    "s",
    "t",
    "u",
    "v",
    "w",
    "x",
    "y",
    "z",
]
SIGNIF_LINES_OFFSET = 0.05


def posthoc_non_parametric(groups, values, alpha):
    # Create a df for the dunn test
    group_col = []
    value_col = []
    for idx in range(len(groups)):
        group_col.extend([groups[idx] for _ in range(len(values[idx]))])
        value_col.extend(values[idx])

    data = pd.DataFrame({'group': group_col, 'value': value_col})

    #print(groups)
    #print(values)

    # Perform Dunn's test with Bonferroni correction
    dunn_result = sp.posthoc_dunn(data, val_col='value', group_col='group', p_adjust='fdr_bh')

    significant_pairs = []

    #print("\nPost-hoc pairwise comparison (Dunn's test):\n")

    for i, j in combinations(groups, 2):
        p = dunn_result[i][j]
        #print(f"Comparison: {i} vs {j}")
        #print(f"P-Value: {p}")
        if p < alpha:
            #print("The difference is statistically significant.\n")
            significant_pairs.append((i, j, p))
        #else:
            #print("The difference is not statistically significant.\n")

    if not significant_pairs:
        significant_pairs = [
            tuple(),
        ]
    
    return significant_pairs


def _shapiro_test(values, alpha):
    passed_normality = []
    p_values = []

    for value in values:
        W, p_value = stats.shapiro(value)
        passed_normality.append(p_value > alpha)
        p_values.append(p_value)

    return p_values, passed_normality


def _bartlett_test(values, alpha):
    W, p_value = stats.bartlett(*values)
    passed_homogeneity = p_value > alpha

    return p_value, passed_homogeneity


def _test_assumptions_and_transform(groups, values, alpha, override_transformation = None):
    transformation_success = False
    
    if override_transformation:
        transformed_values = []
        for data in values:
            transformed_values.append([override_transformation(x) for x in data])
            transformation_success = True
    else:
        normality = _shapiro_test(values, alpha)
        homogeneity = _bartlett_test(values, alpha)
        
        if np.all(normality[1]):
            if homogeneity[1]:
                    print('No transformation needed')
                    return True, values
    
        for idx, transformation in enumerate(TRANSFORMATIONS):
            try:
                transformed_values = []
                for data in values:
                    transformed_values.append([transformation(x) for x in data])
                normality = _shapiro_test(transformed_values, alpha)
                homogeneity = _bartlett_test(transformed_values, alpha)
                if np.all(normality[1]):
                    if homogeneity[1]:
                        transformation_success = True
                        if idx:
                            print("Data transformed with transformation number " + str(idx))
                        break
            except (ValueError, ZeroDivisionError):
                continue

    return transformation_success, transformed_values


def _tukey_test(groups, transformed_values, alpha):
    df_values = []
    df_groups = []
    for idx in range(len(groups)):
        for value in transformed_values[idx]:
            df_values.append(value)
            df_groups.append(groups[idx])
    df = pd.DataFrame({"group": df_groups, "value": df_values})
    tukey_results = pairwise_tukeyhsd(df["value"], df["group"], alpha=alpha)
    return tukey_results, df


def signif_groups_from_tukey(tukey_results):
    tukey_results = tukey_results[0]
    significant_pairs = []
    count = 0
    for result in tukey_results.data[1:]:
        group1, group2, _, p, _, _, reject = result
        if reject:
            significant_pairs.append((group1, group2, p))
            count += 1

    return significant_pairs


def stats_test(groups, values, alpha, transformation = None):
    if not type(values[0]) == np.ndarray:
        values = [np.array(x) for x in values]

    transformation_success, transformed_values = _test_assumptions_and_transform(
        groups, values, alpha, transformation
    )

    if not transformation_success:
        H_statistic, p_value = stats.kruskal(*values)
        #normality = _shapiro_test(values, alpha)
        #homogeneity = _bartlett_test(values, alpha)
        print("Data failed ANOVA assumptions:")
        print(
            "Kruskal Wallis Test Results: H = "
            + str(round(H_statistic, 2))
            + ", p = "
            + str(p_value)
        )
        #print("Groups: " + str(list(groups)))
        #print(
        #    "Normality Shapiro-Wilk test: p = "
        #    + str([round(x, 3) for x in normality[0]])
        #)
        #print("Homogeneity: Levene test: p = " + str(homogeneity[0]))
        #print("Group variances: " + str([round(x.var(), 2) for x in values]) + "\n")
        return (H_statistic, p_value), posthoc_non_parametric(groups, values, alpha)

    f_statistic, p_value = stats.f_oneway(*transformed_values)
    
    print(f'ANOVA results: F = {f_statistic}, p = {p_value}')

    tukey_results_and_df = list(_tukey_test(groups, transformed_values, alpha))
    tukey_results_and_df[0] = tukey_results_and_df[0].summary()
    
    print(tukey_results_and_df[0])

    return (f_statistic, p_value), signif_groups_from_tukey(tukey_results_and_df)
