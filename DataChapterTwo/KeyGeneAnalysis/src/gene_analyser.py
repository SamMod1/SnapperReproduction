import pandas as pd
import sys
import src.stats as stats
import plotly
import plotly.graph_objects as go


SIGNIF_LINES_OFFSET = 0.05
FIG_WIDTH = 500
FIG_HEIGHT = 400


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


def build_signif_symbol(alpha_signif_symbols, p):
    levels = sorted(list(alpha_signif_symbols.keys()), key = lambda x: float(x))
    for signif_level in levels:
        if p <= float(signif_level):
            return alpha_signif_symbols[signif_level]
    
    raise ValueError("Value is not statistically significant. The alpha level or alpha significance levels symbol mappings may be wrong")
    

def annotate_boxplot(fig, unique_groups, values, post_hoc_results, alpha_signif_symbols):
    fig.update_layout(
        showlegend=False,
        width=FIG_WIDTH,
        height=FIG_HEIGHT,
        margin=dict(l=70, r=0, b=40, t=20),
    )
    if type(post_hoc_results[0]) == tuple:
        significant_pairs = post_hoc_results
    else:
        significant_pairs = signif_groups_from_tukey(post_hoc_results)
    if not significant_pairs:
        return fig
    if not significant_pairs[0]:
        return fig

    # Adjust the y_max calculation for plotly
    y_max = max([max(val) for val in values])
    offset = y_max * SIGNIF_LINES_OFFSET

    shapes = []
    annotations = []

    for group1, group2, p in significant_pairs:
        i, j = unique_groups.index(group1), unique_groups.index(group2)
        y = y_max + offset

        # Add lines and asterisk as annotations in plotly
        shapes.append(
            dict(
                type="line",
                xref="x",
                yref="y",
                x0=i,
                y0=y + offset,
                x1=j,
                y1=y + offset,
                line=dict(color="black", width=1.5),
            )
        )
        shapes.append(
            dict(
                type="line",
                xref="x",
                yref="y",
                x0=i,
                y0=y,
                x1=i,
                y1=y + offset,
                line=dict(color="black", width=1.5),
            )
        )
        shapes.append(
            dict(
                type="line",
                xref="x",
                yref="y",
                x0=j,
                y0=y,
                x1=j,
                y1=y + offset,
                line=dict(color="black", width=1.5),
            )
        )
        annotations.append(
            dict(
                x=((i + j + 2) / 2) - 1,
                y=y + offset * 1.25,
                xref="x",
                yref="y",
                text=build_signif_symbol(alpha_signif_symbols, p),
                showarrow=False,
                font=dict(size=12),
            )
        )
        y_max = max(y_max, y + 2 * offset)

    fig.update_layout(annotations=annotations, shapes=shapes)
    return fig


def create_boxplot(x, y, origional_samples, chart_title):
    fig = go.Figure()
    for group, yi, label in zip(x, y, origional_samples):
        if group == 'Female':
            fig.add_trace(go.Box(y=yi, name=group, boxpoints="all", text=label, fillcolor='rgb(246, 169, 156)', line={'color':'rgb(239, 80, 53)'}))
        elif group == 'Male':
            fig.add_trace(go.Box(y=yi, name=group, boxpoints="all", text=label, fillcolor='rgb(176, 182, 252)', line={'color':'rgb(100, 111, 250)'}))
        elif group == 'Juvenile':
            fig.add_trace(go.Box(y=yi, name=group, boxpoints="all", text=label, fillcolor='rgb(192, 192, 192)', line={'color':'rgb(128, 128, 128)'}))
        else:
            #fig.add_trace(go.Box(y=yi, name=group, boxpoints=False, text=label))
            fig.add_trace(go.Box(y=yi, name=group, boxpoints="all", text=label))
    fig.update_layout(yaxis_title="Gene Expression", title=chart_title)
    #fig.update_layout(
    #    paper_bgcolor='rgba(0,0,0,0)',
    #    plot_bgcolor='rgba(0,0,0,0)',
    #)

    return fig


def aggregate_groups(transcript_df):
    groups = sorted(transcript_df["sample"].unique())
    values = [
        [x for x in transcript_df.loc[transcript_df["sample"] == y, "tpm"]]
        for y in groups
    ]
    origional_samples = [
        [x for x in transcript_df.loc[transcript_df["sample"] == y, "origional_sample"]]
        for y in groups
    ]

    return groups, values, origional_samples


def prep_data(df, analysis_levels):
    df["origional_sample"] = pd.Series(df["sample"])
    for level in analysis_levels:
        df.loc[df["sample"].str.contains(level), "sample"] = level
    df = df.loc[df["sample"].isin(analysis_levels)]

    return df


def select_data(gene_df, gene):
    print(f"""{str(len(gene_df))} transcripts for the gene {gene} were found:""")
    for idx in range(len(gene_df)):
        print(str(idx + 1) + ": " + gene_df.iloc[idx]["tx"])
    selection = input(
        """Enter the number(s) of the transcript(s) you want to look at, (seperated by commas)
Alternatively, leave blank to analyse all genes
Or, type 'q' to exit and analyse no genes

Input: """
    )
    if selection == "q":
        sys.exit(0)
    elif len(selection) > 0:
        selection = selection.split(",")
        selection = [int(x) - 1 for x in selection]
        gene_df = gene_df.iloc[selection]

    return gene_df


def drop_empty_cols(df):
    for col in df.columns:
        if df[col].sum() == 0:
            df.drop(col, axis=1, inplace=True)

    return df


def init_data(idx, df):
    gene_df = df.iloc[[idx]]
    gene_df = gene_df.transpose()
    gene_df.columns = gene_df.iloc[0]
    gene_df = gene_df[2:]
    gene_df["sample"] = gene_df.index
    gene_df.index = [x for x in range(len(gene_df))]

    return gene_df


def analyse_gene(gene_df, alpha, gene, plotly_js, alpha_signif_symbols):
    for col_idx in range(len(gene_df.columns) - 2):
        print(gene_df.columns[col_idx] + ":")
        transcript_df = pd.DataFrame(
            {
                "sample": gene_df["sample"],
                "tpm": gene_df[gene_df.columns[col_idx]],
                "origional_sample": gene_df["origional_sample"],
            }
        )
        groups, values, origional_samples = aggregate_groups(transcript_df)
        stats_results = stats.stats_test(groups, values, alpha)
        boxplot = create_boxplot(groups, values, origional_samples, "")
        anova_results, post_hoc_results = stats_results
        if type(post_hoc_results[0]) == tuple:
            print(
                "Kruskal Wallis Results: H = "
                + str(round(anova_results[0], 2))
                + ", p = "
                + str(round(anova_results[1], 2))
            )
            print("")
        else:
            print(
                "Anova Results: f = "
                + str(round(anova_results[0], 2))
                + ", p = "
                + str(round(anova_results[1], 2))
            )
            print(post_hoc_results[0])
            print("")
        boxplot = annotate_boxplot(boxplot, groups, values, post_hoc_results, alpha_signif_symbols)
        boxplot = plotly.offline.plot(
            boxplot, include_plotlyjs=plotly_js, output_type="div"
        )

        return boxplot, anova_results, post_hoc_results


def process_data(
        data, 
        gene, 
        tpm_data, 
        analysis_levels, 
        alpha=0.05, 
        alpha_signif_symbols = {'0.05': '*', '0.01': '**', '0.001': '***'}
    ):
    plotly_js = False#True
            
    for tissue_idx, tissue_df in enumerate(tpm_data):
        gene_df_whole = tissue_df.loc[tissue_df["gene_id"] == gene]
        if len(gene_df_whole) > 1:
            numerical = gene_df_whole[gene_df_whole.columns[2:]]
            sums = numerical.sum()
            gene_df_whole.loc[-1] = list([gene, gene] + list(sums))
            gene_df_whole.index = gene_df_whole.index + 1  # shifting index
            gene_df_whole.sort_index(inplace=True)
        for transcript_idx in range(len(gene_df_whole)):
            gene_df = init_data(transcript_idx, gene_df_whole)
            samples_saved = list(gene_df["sample"])
            transcript_id = gene_df.columns[0]
            data[transcript_idx].append(transcript_id)
            for analysis_level_idx, analysis_level in enumerate(analysis_levels):
                analysis_df = prep_data(pd.DataFrame(gene_df), analysis_level)
                gene_df["sample"] = samples_saved
                analysis_df = drop_empty_cols(analysis_df)
                boxplot = analyse_gene(analysis_df, alpha, gene, plotly_js, alpha_signif_symbols)
                # if boxplot:
                #     plotly_js = False
                data[transcript_idx][tissue_idx][analysis_level_idx] = boxplot

    return data


def find_gene(gene, tpm_data, analysis_levels, exact_match = True):
    if exact_match:
        n_transcripts = max([len(x.loc[x["gene_id"] == gene]) for x in tpm_data])
        unique_genes = [gene]
    else:
        genes = [x.loc[x["gene_id"].str.startswith(gene)]['gene_id'] for x in tpm_data]
        n_transcripts = max([len(x) for x in genes])
        unique_genes = [list(x.unique()) for x in genes]
        unique_genes = set([
            x
            for xs in genes
            for x in xs
        ])
    
    if not n_transcripts:
        return None, None
    
    if n_transcripts == 1:
        data = [[[["NONE"] for _ in range(len(analysis_levels))] for _ in range(len(tpm_data))]]
    else:
        data = [
            [[["NONE"] for _ in range(len(analysis_levels))] for _ in range(len(tpm_data))]
            for _ in range(n_transcripts + 1)
        ]
    
    return data, unique_genes

def gene_selection(genes):
    genes = list(genes)
    while True:
        text = 'Multiple genes were found. Please select the one to analyse or enter nothing to analyse all:\n'
        
        for count, gene in enumerate(genes):
            text += str(count + 1) + ': ' + gene + '\n'
        
        text += '\nEnter gene number: '
        gene = input(text)
        
        if gene == 'q':
            return
        elif not gene:
            return genes
            
        try:
            gene = genes[int(gene) - 1]
            return [gene]
        except ValueError:
            print("Error: Entry was not a number, enter 'q' if you wish to quit\n\n")


def build_html(
    gene_name, data, main_template, tab_template, button_str, plotly_js, tpm_files, alpha
):
    main_template = main_template.replace("$gene_name", gene_name)

    for transcript_idx, transcript in enumerate(data):
        this_tab = str(tab_template)
        for tissue_idx, tissue in enumerate([x for x in transcript if type(x) == list]):
            this_tab = this_tab.replace(
                f"$test_info_{str(tissue_idx)}", ""
            )
            for al_idx, analysis_level in enumerate(tissue):
                try:
                    this_tab = this_tab.replace(
                        f"$al_{str(al_idx)}_tissue_{str(tissue_idx)}",
                        data[transcript_idx][tissue_idx][al_idx][0],
                    )
                except TypeError:
                    this_tab = this_tab.replace(
                        f"$al_{str(al_idx)}_tissue_{str(tissue_idx)}", ""
                    )
        this_tab = this_tab.replace("$transcript_id", str(transcript_idx))
        this_tab = this_tab.replace("$transcript", data[transcript_idx][-1])
        main_template = main_template.replace("$new_tab", this_tab)
        button = button_str.replace("$transcript_id", str(transcript_idx))
        button = button.replace("$transcript", data[transcript_idx][-1])
        main_template = main_template.replace("$new_transcript_button", button)

    main_template = main_template.replace("$new_tab", "")
    main_template = main_template.replace(
        "$data_source", ", ".join([x for x in tpm_files])
    )
    main_template = main_template.replace("$alpha_value", str(alpha))
    main_template = main_template.replace("$new_transcript_button", "")
    main_template = main_template.replace('$plotlyjs', plotly_js)

    return main_template


def run_analysis(
    gene,
    tpm_files,
    analysis_levels,
    main_template,
    tab_template,
    button_string,
    plotly_js,
    alpha,
    alpha_signif_symbols,
    data,
    tpm_data
):
    data = process_data(data, gene, tpm_data, analysis_levels, alpha, alpha_signif_symbols)

    gene = gene.replace('gene-', '')
    html = build_html(
        gene, data, main_template, tab_template, button_string, plotly_js, tpm_files, alpha
    )
    gene = gene.replace(':', ';')
    open(f"output//{gene}.html", "w", encoding="utf-8").write(html)


def main(
    tpm_files,
    analysis_levels,
    main_template,
    tab_template,
    button_string,
    plotly_js,
    samples_to_exclude = None,
    alpha = 0.05,
    alpha_signif_symbols = {'0.05': '*', '0.01': '**', '0.001': '***'},
    labels = None
):
    if not labels:
        labels = analysis_levels
    plotly_js = open(plotly_js, 'r', encoding = 'utf-8').read()
    main_template = open(main_template, "r", encoding="utf-8").read()
    tab_template = open(tab_template, "r", encoding="utf-8").read()
    tpm_data = [pd.read_csv(x, sep="\t") for x in tpm_files]
    
    if samples_to_exclude:
        for idx in range(len(tpm_data)):
            tpm_data[idx] = tpm_data[idx][[x for x in tpm_data[idx].columns if not x in samples_to_exclude]]
    
    for df in tpm_data:
        new_columns = list(df.columns)
        for idx_1, level_1 in enumerate(analysis_levels):
            for idx_2, level_2 in enumerate(level_1):
                for col_idx in range(len(new_columns)):
                    new_columns[col_idx] = new_columns[col_idx].replace(
                        level_2, labels[idx_1][idx_2]
                    )
        df.columns = new_columns
    analysis_levels = labels
        
    while True:
        gene = input("Name of gene to analyse (or 'q' to quit): ")
        if gene == "q":
            return
        if not gene.startswith('gene-'):
            gene = 'gene-' + gene
            
        data, _ = find_gene(gene, tpm_data, analysis_levels)
        if not data:
            data, unique_genes = find_gene(gene, tpm_data, analysis_levels, exact_match = False)
            if data:
                genes = gene_selection(unique_genes)
                if not genes:
                    continue
                for gene in genes:
                    data, _ = find_gene(gene, tpm_data, analysis_levels)
                    run_analysis(
                        gene,
                        tpm_files,
                        analysis_levels,
                        main_template,
                        tab_template,
                        button_string,
                        plotly_js,
                        alpha,
                        alpha_signif_symbols,
                        data,
                        tpm_data
                    )
                    
            if not data:
                print('No genes found for: ' + gene)
                continue
        
        run_analysis(
            gene,
            tpm_files,
            analysis_levels,
            main_template,
            tab_template,
            button_string,
            plotly_js,
            alpha,
            alpha_signif_symbols,
            data,
            tpm_data
        )
