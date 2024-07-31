import sys
import ast
import random
import pandas as pd
import skbio
import seaborn as sns
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import scikit_posthocs as sp
import skbio.diversity.alpha as alpha
from skbio import TreeNode, DistanceMatrix
from skbio.tree import nj
from skbio.diversity import beta_diversity
from skbio.stats.ordination import pcoa
from skbio.stats.distance import permanova
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.manifold import MDS
from matplotlib import cm
from matplotlib.patches import Ellipse
from matplotlib_venn import venn2, venn3
from scipy.stats import f_oneway, shapiro, mannwhitneyu, kruskal, spearmanr, pearsonr
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from sklearn import preprocessing
from pydeseq2.preprocessing import deseq2_norm


def one_condition_struct(asv_info, condition): # Function to create a structure containing different conditions with associated samples
    # Check if a condition file is provided
    if condition is None:
        print("Error: Please provide a condition file.")
        return

    print("Name of conditions available :")
    available_condition = list(condition.columns)
    for name in available_condition:
        print(name)
    column_sample_name = input("What is the name of the column with the sample names? : ")
    column_condition_name = input("What is the name of the column with the names of the conditions to which the samples belong? : ")
    
    # Create a dictionary to store the conditions and their associated samples
    conditions = {}
    for sample, cond in zip(condition[column_sample_name], condition[column_condition_name]):
        matching_samples = [col for col in asv_info.columns if col.startswith(sample)] # Find columns in asv_info DataFrame that match the sample name
        
        if not matching_samples:
            print(f"No matching samples found for {sample}. Skipping.")
            continue
            
        # If the condition is not already in the dictionary, add it with an empty list : allows conditions to be written only once
        if cond not in conditions:
            conditions[cond] = []
        
        # Extend the list of samples associated with the condition
        conditions[cond].extend(matching_samples)

    return conditions


def two_conditions_struct(asv_info, condition): # Function to create a structure containing different conditions with associated samples
    # Check if a condition file is provided
    if condition is None:
        print("Error: Please provide a condition file.")
        return

    print("Name of conditions available :")
    available_condition = list(condition.columns)
    for name in available_condition:
        print(name)
    column_sample_name = input("What is the name of the column with the sample names? : ")
    column_condition_first_name = input("What is the name of the column with the names of the first condition to which the samples belong? :  ")
    column_condition_second_name = input("What is the name of the column with the names of the second condition to which the samples belong? :  ")
    
    # Create a dictionary to store the conditions and their associated samples
    conditions = {}
    for sample, cond1, cond2 in zip(condition[column_sample_name], condition[column_condition_first_name], condition[column_condition_second_name]):
        matching_samples = [col for col in asv_info.columns if col.startswith(sample)] # Find columns in asv_info DataFrame that match the sample name
        
        if not matching_samples:
            print(f"No matching samples found for {sample}. Skipping.")
            continue
            
        # If conditions are not already in the dictionary, add them with an empty list : allows conditions to be written only once
        if (cond1, cond2) not in conditions:
            conditions[(cond1, cond2)] = []
            
        # Extend the list of samples associated with the condition
        conditions[(cond1, cond2)].extend(matching_samples)
        
    print(conditions)

    return conditions, column_condition_first_name, column_condition_second_name


#####################
## Alpha Diversity ##
#####################


def alpha_diversity_one(asv_info): # Function to calculate alpha diversity for a given sample
    print("Available samples :")
    for column in asv_info.columns[1:]:
        print(column)
    sample_alpha = input("Which sample do you want alpha diversity ? : ")
    counts = asv_info[sample_alpha]     # Retrieve counts for the selected sample
    alpha_index = input("Which alpha diversity index do you want to calculate ? (shannon / simpson / inverse_simpson / chao / richness): ")
    if alpha_index == 'shannon':
        alpha_diversity = skbio.diversity.alpha.shannon(counts)
        print("-- Alpha Diversity : Shannon index for sample ", sample_alpha, " : ", alpha_diversity)
        
    elif alpha_index == 'simpson':
        alpha_diversity = skbio.diversity.alpha.simpson(counts)
        print("-- Alpha Diversity : Simpson index for sample ", sample_alpha, " : ", alpha_diversity)
        
    elif alpha_index == 'inverse_simpson':
        simpson_index = skbio.diversity.alpha.simpson(counts)
        alpha_diversity = 1 / simpson_index if simpson_index != 0 else float('inf') # Calculate the inverse Simpson index
        print("-- Alpha Diversity : Inverse Simpson index for sample ", sample_alpha, " : ", alpha_diversity)
            
    elif alpha_index == 'chao':
        alpha_diversity = skbio.diversity.alpha.chao1(counts)
        print("-- Alpha Diversity : Chao index for sample ", sample_alpha, " : ", alpha_diversity)
        
    elif alpha_index == 'richness':
        asv_sample = asv_info[['ASVNumber', sample_alpha]].loc[asv_info[sample_alpha] > 0]
        asv_names = asv_sample['ASVNumber'].tolist()
        asv_number = len(asv_names)
        print("-- Alpha Diversity : Observed richness for sample ",sample_alpha, " : ",  asv_number)
    
    else:
        print("Alpha diversity index not supported.")
        exit()
    
   
def alpha_diversity_all(asv_info): # Function to calculate alpha diversity for all samples
    alpha_diversity_all = {}
    alpha_index = input("Which alpha diversity index do you want to calculate ? (shannon / simpson / inverse_simpson / chao / richness): ")
    if alpha_index == 'shannon':
        print("-- Shannon Alpha diversity for all samples --")
        for column in asv_info.columns[1:]:  # Ignore the first column (ASVNumber) to have samples columns
            counts = asv_info[column] # Retrieve counts for the selected sample
            alpha_diversity_all[column] = {}
            alpha_diversity_all[column]['shannon'] = skbio.diversity.alpha.shannon(counts)
        for column, diversity in alpha_diversity_all.items():
            print(column, " : ", diversity['shannon'])
            
    elif alpha_index == 'simpson':
        print("-- Simpson Alpha diversity for all samples --")
        for column in asv_info.columns[1:]:
            counts = asv_info[column]
            alpha_diversity_all[column] = {}
            alpha_diversity_all[column]['simpson'] = skbio.diversity.alpha.simpson(counts)
        for column, diversity in alpha_diversity_all.items():
            print(column, " : ", diversity['simpson'])
            
    elif alpha_index == 'inverse_simpson':
        print("-- Inverse Simpson Alpha diversity for all samples --")
        for column in asv_info.columns[1:]:
            counts = asv_info[column]
            alpha_diversity_all[column] = {}
            simpson_index = skbio.diversity.alpha.simpson(counts)
            inverse_simpson_index = 1 / simpson_index if simpson_index != 0 else float('inf') # Calculate the inverse Simpson index
            alpha_diversity_all[column]['inverse_simpson'] = inverse_simpson_index
        for column, diversity in alpha_diversity_all.items():
            print(column, " : ", diversity['inverse_simpson'])
            
    elif alpha_index == 'chao':
        print("-- Chao Alpha diversity for all samples --")
        for column in asv_info.columns[1:]:
            counts = asv_info[column]
            alpha_diversity_all[column] = {}
            alpha_diversity_all[column]['chao'] = skbio.diversity.alpha.chao1(counts)
        for column, diversity in alpha_diversity_all.items():
            print(column, " : ", diversity['chao'])
            
    elif alpha_index == 'richness':
        print("-- Observed Richness Alpha diversity for all samples --")
        asv_all_samples = {}
        for column in asv_info.columns[1:]:
            asv_samples = asv_info[['ASVNumber', column]].loc[asv_info[column] > 0] # Keep only rows of the df where the value in the current column is greater than 0
            asv_names = asv_samples['ASVNumber'].tolist()
            asv_all_samples[column] = len(asv_names) # Number of richness ASVs for the current sample
        for sample, asv in asv_all_samples.items():
            print(sample, " : ", asv)
    else:
        print("Alpha diversity index not supported.")
        exit()
    

def statistical_test_alpha(asv_info, condition): # Function to perform statistical tests on alpha diversity
    # Retrieve the structure containing sample data based on conditions
    conditions=one_condition_struct(asv_info, condition)
    
    # Calculate alpha diversity for each condition
    alpha_results = {}
    alpha_index = input("On which alpha diversity index would you like to perform a statistical test ? (shannon / simpson / inverse_simpson / chao / richness): ")

    for cond, samples in conditions.items():
        alpha_results[cond] = []
        for sample in samples:
            counts = asv_info[sample]
            if alpha_index == 'shannon':
                alpha_results[cond].append(skbio.diversity.alpha.shannon(counts))
            elif alpha_index == 'simpson':
                alpha_results[cond].append(skbio.diversity.alpha.simpson(counts))
            elif alpha_index == 'inverse_simpson':
                simpson_index = skbio.diversity.alpha.simpson(counts)
                alpha_results[cond].append(1 / simpson_index if simpson_index != 0 else float('inf'))
            elif alpha_index == 'chao':
                alpha_results[cond].append(skbio.diversity.alpha.chao1(counts))
            elif alpha_index == 'richness':
                asv_samples = asv_info[['ASVNumber', sample]].loc[asv_info[sample] > 0]
                asv_names = asv_samples['ASVNumber'].tolist()
                alpha_results[cond].append(len(asv_names))
                
    print(alpha_results)
    
    alpha_values = []
    for values in alpha_results.values(): # Boucle sur les valeurs de alpha_results
        alpha_values.extend(values)
          
    s_statistic, p_value = shapiro(alpha_values)
    print("Shapiro-Wilk test result for normality (alpha diversity) : ")
    print("Statistic for the test : ", s_statistic)
    print("p-value : ", p_value)
    
    threshold_test = input("Threshold for test : ")
    
    if p_value > float(threshold_test):
        print(" Alpha diversity values comes from a normal distribution  : ")
        test_para = input("Do you want to make ANOVA test or Pearson correlation test ? (anova/pearson) : ")
        
        if test_para == 'anova':
            alpha_values_array = [np.array(values) for values in alpha_results.values()]
            print(alpha_values_array)
            # ANOVA test
            f_statistic, p_value = f_oneway(*alpha_values_array) # One-way ANOVA test
            print("ANOVA test result between all groups : ")
            print("Statistic F : ", f_statistic)
            print("p-value : ", p_value)
                    
            post_hoc = input("Do you want to make post-hoc test (Tukey) ? Yes/No : ") # 95% test
            if post_hoc == 'Yes':
                alpha_values_flat = np.concatenate(alpha_values_array) # Flattening the data
                conditions_list = [] # Creating a list of conditions for Tukey test
                for cond, data in alpha_results.items():
                    for value in data:
                        conditions_list.append(cond) # Add conditions for each sample

                tukey_result = pairwise_tukeyhsd(endog=alpha_values_flat, groups=conditions_list)
                print("Result of Tukey test : ")
                print(tukey_result)
            
        elif test_para == 'pearson':
            # Pearman correlation test
            print("Name of groups available :")
            for cond, samples in conditions.items():
                print(cond)

            first_group = input("Name of first group you want to compare samples : ")
            second_group = input("Name of second group you want to compare samples : ")
                        
            first_array = [values for values in alpha_results[first_group]]
            second_array = [values for values in alpha_results[second_group]]
            
            if len(first_array) != len(second_array): # All the input array must have the same length.
                min_length = min(len(first_array), len(second_array)) # Crop the array to the min length array
                first_array = first_array[:min_length]
                second_array = second_array[:min_length]
            
            corr_statistic, p_value = pearsonr(first_array, second_array)
            print("Pearson correlation test result between these two groups : ")
            print("Correlation coefficient : ", corr_statistic)
            print("p-value : ", p_value)
        
        else:
            print("Statistical test not supported.")
        
                
    else:
        print(" Alpha diversity values does not come from a normal distribution  : ")
        test_non_para = input("Do you want to make Kruskal-Wallis test or Mann-Whitney test or Spearman correlation test ? (kruskal / mann / spearman) : ")
        
        if test_non_para == 'kruskal':
            alpha_values_array = [np.array(values) for values in alpha_results.values()]
            # Kruskal-Wallis test
            h_statistic, p_value = kruskal(*alpha_values_array)
            print("Kruskal-Wallis test result between all groups : ")
            print("Statistic H : ", h_statistic)
            print("p-value : ", p_value)
                    
            post_hoc = input("Do you want to make post-hoc test (Dunn) ? Yes/No : ")
            if post_hoc == 'Yes':
                alpha_values_flat = np.concatenate(alpha_values_array) # Flattening the data
            
                conditions_list = [] # Creating a list of conditions for Tukey test
                for cond, data in alpha_results.items():
                    for value in data:
                        conditions_list.append(cond) # Add conditions for each sample
                
                df = pd.DataFrame({'Values': alpha_values_flat,'Groups': conditions_list})
                
                dunn_result = sp.posthoc_dunn(df, val_col='Values', group_col='Groups')
                print("Result of Dunn test : ")
                print(dunn_result)
    
        elif test_non_para == 'mann':
            print("Name of groups available :")
            for cond, samples in conditions.items():
                print(cond)

            first_group = input("Name of first group you want to compare samples : ")
            second_group = input("Name of second group you want to compare samples : ")
                        
            first_array = [values for values in alpha_results[first_group]]
            second_array = [values for values in alpha_results[second_group]]
            
            mwu_statistic, p_value = mannwhitneyu(first_array, second_array) # Mann-Whitney U test
            print("Mann-Whitney U test result between these two groups : ")
            print("Statistic for the test : ", mwu_statistic)
            print("p-value : ", p_value)
            
        elif test_non_para == 'spearman':
            # Spearman correlation test
            print("Name of groups available :")
            for cond, samples in conditions.items():
                print(cond)

            first_group = input("Name of first group you want to compare samples : ")
            second_group = input("Name of second group you want to compare samples : ")
                        
            first_array = [values for values in alpha_results[first_group]]
            second_array = [values for values in alpha_results[second_group]]
        
            if len(first_array) != len(second_array): # All the input array dimensions must match exactly
                min_length = min(len(first_array), len(second_array)) # Crop the array to the min length array
                first_array = first_array[:min_length]
                second_array = second_array[:min_length]
            
            corr_statistic, p_value = spearmanr(first_array, second_array)
            print("Spearman correlation test result between these two groups : ")
            print("Correlation coefficient : ", corr_statistic)
            print("p-value : ", p_value)
        
        else:
            print("Statistical test not supported.")
        

def alpha_graph(asv_info): # Function to create alpha diversity scatter plot
    alpha_index = input("Which alpha diversity index would you like ? (shannon / simpson / inverse_simpson / chao / richness) : ")
    alpha_diversity = {} # Dictionary to store alpha diversity values for each sample and index
    for column in asv_info.columns[1:]:
        counts = asv_info[column] # Extract counts for the sample
        if alpha_index == 'shannon':
            alpha_diversity.setdefault('shannon', {})[column] = skbio.diversity.alpha.shannon(counts)
        elif alpha_index == 'simpson':
            alpha_diversity.setdefault('simpson', {})[column] = skbio.diversity.alpha.simpson(counts)
        elif alpha_index == 'inverse_simpson':
            simpson_index = skbio.diversity.alpha.simpson(counts)
            alpha_diversity.setdefault('inverse_simpson', {})[column] = 1 / simpson_index if simpson_index != 0 else float('inf')
        elif alpha_index == 'chao':
            alpha_diversity.setdefault('chao', {})[column] = skbio.diversity.alpha.chao1(counts)
        elif alpha_index == 'richness':
            asv_samples = asv_info[['ASVNumber', column]].loc[asv_info[column] > 0] # Keep only rows of the df where the value in the current column is greater than 0
            asv_names = asv_samples['ASVNumber'].tolist()
            alpha_diversity.setdefault('richness', {})[column] = len(asv_names) # Number of richness ASVs for the current sample
        else:
            print("Alpha diversity index not supported.")
            exit()

    # Create the figure
    plt.figure(figsize=(12, 10))
    plt.scatter(alpha_diversity[alpha_index].keys(), alpha_diversity[alpha_index].values(), color='mediumturquoise')
    plt.xlabel('Samples')
    plt.ylabel(f'{alpha_index.capitalize()} Alpha Diversity')
    plt.title(f'{alpha_index.capitalize()} Alpha Diversity according to samples')
    plt.xticks(rotation=90)
    plt.tight_layout()
    format_file = input("Which file format would you like to save the plot ? SVG / PDF / PNG ")
    if format_file == 'PDF':
        plt.savefig("scatter_plot_alpha_div.pdf", format='pdf', pad_inches=0.2)
    elif format_file == 'SVG':
        plt.savefig("scatter_plot_alpha_div.svg", format='svg', pad_inches=0.2)
    elif format_file == 'PNG':
        plt.savefig("scatter_plot_alpha_div.png", format='png', pad_inches=0.2)
    plt.show()


def alpha_graph_condition(asv_info, condition): # Function to create alpha diversity boxplots grouped by conditions
    alpha_index = input("Which alpha diversity index would you like? (shannon / simpson / inverse_simpson  / chao / richness): ")
    # Retrieve the structure containing sample data based on conditions
    conditions=one_condition_struct(asv_info, condition)
    
    # Calculate alpha diversity for each condition
    alpha_results = {}
    for cond, samples in conditions.items():
        alpha_results[cond] = []
        for sample in samples:
            counts = asv_info[sample]
            if alpha_index == 'shannon':
                alpha_results[cond].append(skbio.diversity.alpha.shannon(counts))
            elif alpha_index == 'simpson':
                alpha_results[cond].append(skbio.diversity.alpha.simpson(counts))
            elif alpha_index == 'inverse_simpson':
                simpson_index = skbio.diversity.alpha.simpson(counts)
                alpha_results[cond].append(1 / simpson_index if simpson_index != 0 else float('inf'))
            elif alpha_index == 'chao':
                alpha_results[cond].append(skbio.diversity.alpha.chao1(counts))
            elif alpha_index == 'richness':
                asv_samples = asv_info[['ASVNumber', sample]].loc[asv_info[sample] > 0]
                asv_names = asv_samples['ASVNumber'].tolist()
                alpha_results[cond].append(len(asv_names))
            else:
                print("Alpha diversity index not supported.")
                exit()

    # Set the color palette
    colors = sns.color_palette('rainbow', n_colors=len(alpha_results))
    
    # Boxplots for alpha diversity grouped by conditions
    plt.figure(figsize=(10, 6))
    bp = plt.boxplot(alpha_results.values(), labels=alpha_results.keys(), patch_artist=True, capprops={'linewidth': 0.0})
    plt.xlabel('Conditions')
    plt.ylabel(f'{alpha_index.capitalize()} Alpha Diversity')
    plt.title(f'{alpha_index.capitalize()} Alpha Diversity by Condition')
    plt.xticks(rotation=90, ha='right')
        
    # Customize boxplot and whisker colors and median line
    for i, (box, median, color) in enumerate(zip(bp['boxes'], bp['medians'], colors)):
        box.set(facecolor=color + (0.2,), edgecolor=color)  # Set box color
        median.set(color=color)  # Set median line color
        bp['whiskers'][i * 2].set(color=color)  # Lower whisker
        bp['whiskers'][(i * 2) + 1].set(color=color)  # Upper whisker
    
    # Add points for each value
    for i, cond in enumerate(alpha_results.keys()):
        y = alpha_results[cond]
        x = [i + 1] * len(y)
        plt.plot(x, y, 'k.', alpha=0.9, markersize=9, color=colors[i])
    plt.xticks(ticks=np.arange(1, len(alpha_results) + 1), labels=alpha_results.keys(), rotation=90, ha='right')
    plt.tight_layout()
    format_file = input("Which file format would you like to save the plot ? SVG / PDF / PNG ")
    if format_file == 'PDF':
        plt.savefig("boxplot_alpha_div.pdf", format='pdf', pad_inches=0.2)
    elif format_file == 'SVG':
        plt.savefig("boxplot_alpha_div.svg", format='svg', pad_inches=0.2)
    elif format_file == 'PNG':
        plt.savefig("boxplot_alpha_div.png", format='png', pad_inches=0.2)
    plt.show()



####################
## Beta Diversity ##
####################


def beta_diversity_function(asv_info): # Function to calculate beta diversity
    asv_info.index = asv_info['ASVNumber']
    asv_info = asv_info.drop(columns=['ASVNumber'])
    asv_info = asv_info.T.fillna(0) # Transpose so that samples are rows
        
    beta_index = input("Which alpha diversity indices would you like ? : (braycurtis / jaccard / weighted_unifrac / unweighted_unifrac)")
    if beta_index in ['braycurtis', 'jaccard']:
        # Calculate beta diversity using the selected index and count table
        beta_diversity = skbio.diversity.beta_diversity(beta_index, asv_info)
    
    elif beta_index in ['weighted_unifrac', 'unweighted_unifrac']:

        # Load phylogenetic tree
        tree = TreeNode.read('./phyloTree.newick')
        
        midpoint_rooted_tree = tree.root_at_midpoint() # Set root to the tree
        
        # Function to replace full name with short identifier
        def update_node_names(node):
            if node.name:
                node.name = node.name.split('|')[0] # Extract short identifier (before first '|')
            for child in node.children: # Call recursively on children
                update_node_names(child)

        # Update node names
        update_node_names(midpoint_rooted_tree)
        
        tree_tip_names = {tip.name for tip in midpoint_rooted_tree.tips()} # Get the names of the tips (leaves) of the tree
        otu_ids = set(asv_info.columns) # Obtain the correspondence between the names of the ASVs in the table and the names of the leaves of the tree

        # Check matches
        missing_in_tree = otu_ids - tree_tip_names
        missing_in_otu_table = tree_tip_names - otu_ids
        print(f"OTUs present in the table but missing in the tree : {missing_in_tree}")
        print(f"Tips present in the tree but missing in the table : {missing_in_otu_table}")

        # Calculate UniFrac distances
        beta_diversity = skbio.diversity.beta_diversity(
            metric=beta_index,
            counts=asv_info.values,
            ids=asv_info.index,
            otu_ids=asv_info.columns,
            tree=midpoint_rooted_tree
        )
        
    return beta_diversity
    
    
def beta_diversity_all(beta_diversity): # Function to display beta diversity for all samples
    print("-- Beta diversity distance matrix :")
    print(beta_diversity)

        
def statistical_test_beta(beta_diversity, asv_info, condition): # Function to perform a statistical test (permanova) on beta diversity
    # Retrieve the structure containing sample data based on conditions
    conditions=one_condition_struct(asv_info, condition)
    
    conditions_list = []
    # Extract groups from sample
    for cond, samples in conditions.items():
        for sample in samples:
            conditions_list.append(cond) # Add conditions for each sample
            
    test_diff = input("Do you want to make permanova test on all conditions or Pairwise permanova tests ? (permanova / pairwise) : ")
        
    if test_diff == 'permanova':
        permanova_result = permanova(distance_matrix=beta_diversity, grouping=conditions_list)
        print("-- Permanova test results")
        print("Permanova statistic : ", permanova_result['test statistic'])
        print("p-value : ", permanova_result['p-value'])
        print("Number of permutations : ", permanova_result['number of permutations'])
    
    elif test_diff == 'pairwise':
        unique_cond = set(conditions_list)
        pairwise_results = []

        for cond1 in unique_cond:
            for cond2 in unique_cond:
                if cond1 < cond2:
                    # Filter samples which belong to cond1 et cond2
                    indices = [i for i, cond in enumerate(conditions_list) if cond in [cond1, cond2]]
                    ids = [beta_diversity.ids[i] for i in indices]
                    sub_matrix = beta_diversity.filter(ids)
                    sub_groups = [conditions_list[i] for i in indices]
                    
                    # Pairwise permanova test
                    result = permanova(distance_matrix=sub_matrix, grouping=sub_groups)
                    pairwise_results.append({
                        'Group1': cond1,
                        'Group2': cond2,
                        'Permanova_statistic': result['test statistic'],
                        'p-value': result['p-value'],
                        'Number_of_permutations': result['number of permutations']
                    })

        # Convert into dataframe
        pairwise_results_df = pd.DataFrame(pairwise_results)
        print("-- Pairwise Permanova test results")
        print(pairwise_results_df)
    else:
        print("Statistical test not supported.")
    
    
def beta_diversity_graph(beta_diversity, condition, asv_info): # Function to create beta diversity visualizations
    beta_representation = input("Which beta diversity representation would you like ? (heatmap / NMDS / NMDS_grouped_by_condition / PCoA / PCoA_grouped_by_condition): ")
    
    if beta_representation == 'heatmap': # Create a heatmap representation of beta diversity
        plt.figure(figsize=(10, 8))
        sns.heatmap(beta_diversity.to_data_frame(), cmap="rainbow_r", annot=True, fmt=".2f", linewidths=.5)
        plt.title("Beta diversity heatmap")
        plt.xlabel("Samples")
        plt.ylabel("Samples")
        format_file = input("Which file format would you like to save the plot ? SVG or PDF or PNG ")
        if format_file == 'PDF':
            plt.savefig("heatmap_beta_diversity.pdf", format='pdf', pad_inches=0.2)
        elif format_file == 'SVG':
            plt.savefig("heatmap_beta_diversity.svg", format='svg', pad_inches=0.2)
        elif format_file == 'PNG':
            plt.savefig("heatmap_beta_diversity.png", format='png', pad_inches=0.2)
        plt.show()
        
    elif beta_representation == 'NMDS':
        mds = MDS(metric=False, random_state=0) # MDS object which perform nonmetric MDS with reproducibility
        beta_diversity_array = np.array(beta_diversity[:][:])
        mds_results = mds.fit_transform(beta_diversity_array) # Performs MDS transformation
        stress = mds.stress_
        print(stress)
         
        if 0.1 <= stress <= 0.2:
            mds = MDS(n_components=3, metric=False, random_state=0) # MDS with 3 dimensions
            mds_results = mds.fit_transform(beta_diversity_array) # Performs MDS transformation
            stress = mds.stress_
            print(stress)
            
        if stress <= 0.2:
            colors = sns.color_palette('rainbow', n_colors=len(beta_diversity_array)) # Generates a color for each sample
            mds_results_df = pd.DataFrame(mds_results, columns=['Dimension 1', 'Dimension 2']) # Converts the MDS results into a DataFrame
            sample_names = list(asv_info.columns[1:])
            mds_results_df.index = sample_names # Set the index of the DataFrame to the sample names

            # Create the figure with a legend
            legend_handles = []
            plt.figure(figsize=(12, 8))
            for sample, color in zip(sample_names, colors): # Iterate over sample names and corresponding colors
                # Retrieve the x and y coordinates for the current sample from a DataFrame
                x = mds_results_df.loc[sample, 'Dimension 1']
                y = mds_results_df.loc[sample, 'Dimension 2']
                handle = plt.scatter(x, y, color=color, label=sample) # Plot the current sample with specified color and label
                legend_handles.append(handle) # Append the handle to the legend_handles list for creating the legend later
            plt.legend(handles=legend_handles, title='Samples', loc='best', bbox_to_anchor=(1, 1))
            plt.tight_layout()
            plt.xlabel('NMDS 1')
            plt.ylabel('NMDS 2')
            plt.title('NMDS Plot')
            plt.annotate(f'Stress: {stress:.4f}', xy=(0.83, -0.06), xycoords='axes fraction', fontsize=12,
                 verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.5)) # Add an annotation for the stress value

            format_file = input("Which file format would you like to save the plot ? SVG or PDF or PNG ")
            if format_file == 'PDF':
                plt.savefig("NMDS.pdf", format='pdf', bbox_inches='tight', pad_inches=0.2)
            elif format_file == 'SVG':
                plt.savefig("NMDS.svg", format='svg', bbox_inches='tight', pad_inches=0.2)
            elif format_file == 'PNG':
                plt.savefig("NMDS.png", format='png', bbox_inches='tight', pad_inches=0.2)
            plt.show()
        else :
            print("The stress variable is ", stress, ". It's greater than 0.2. Perform a PCoA analysis instead.")
            pcoa_results = pcoa(beta_diversity) # Perform PCoA
            sample_names = list(asv_info.columns[1:]) # Extract sample name
            pcoa_results.samples.index = sample_names # Set the IDs of PCoA results samples to the sample names
            colors = sns.color_palette('rainbow', n_colors=len(sample_names)) # Generate a list of colors
            legend_handles = []
            plt.figure(figsize=(12, 8))
            
            for sample, color in zip(sample_names, colors): # Iterate over sample names and corresponding colors
                x = pcoa_results.samples.loc[sample, 'PC1']
                y = pcoa_results.samples.loc[sample, 'PC2']
                handle = plt.scatter(x, y, color=color, label=sample)
                legend_handles.append(handle)
                
            plt.legend(handles=legend_handles, title='Samples', loc='best', bbox_to_anchor=(1, 1))
            plt.tight_layout()
            plt.title('PCoA Plot')
            plt.xlabel('Axis 1')
            plt.ylabel('Axis 2')
            format_file = input("Which file format would you like to save the plot ? SVG or PDF or PNG")
            if format_file == 'PDF':
                plt.savefig("PCoA.pdf", format='pdf', bbox_inches='tight', pad_inches=0.2)
            elif format_file == 'SVG':
                plt.savefig("PCoA.svg", format='svg', bbox_inches='tight', pad_inches=0.2)
            elif format_file == 'PNG':
                plt.savefig("PCoA.png", format='png', bbox_inches='tight', pad_inches=0.2)
            plt.show()
        
    elif beta_representation == 'NMDS_grouped_by_condition':
        mds = MDS(metric=False, random_state=0) # MDS object which perform nonmetric MDS with reproducibility
        beta_diversity_array = np.array(beta_diversity[:][:])
        mds_results = mds.fit_transform(beta_diversity_array) # Performs MDS transformation
        stress = mds.stress_
        print(stress)
        
        if 0.1 <= stress <= 0.2:
            mds = MDS(n_components=3, metric=False, random_state=0) # NMDS calculate with 3 dimensions
            mds_results = mds.fit_transform(beta_diversity_array) # Performs NMDS transformation
            stress = mds.stress_
            print(stress)
        
        number_condition = input("How many conditions do you have? (one / two): ")
        if number_condition == 'one':
            if stress <= 0.2:
                mds_results_df = pd.DataFrame(mds_results, columns=['Dimension 1', 'Dimension 2'])
                # Replace row and column numbers with sample names
                sample_names = list(asv_info.columns[1:])
                mds_results_df.index = sample_names
                
                # Retrieve the structure containing sample data based on conditions
                conditions=one_condition_struct(asv_info, condition)
                # Extract legend labels from conditions
                legend_labels = []
                for condition,samples_list in conditions.items():
                    legend_labels.append(condition)
                    
                # Make unique colors from seaborn palette rainbow for each condition
                color = sns.color_palette('rainbow', n_colors=len(legend_labels)) # Generate as many colors as conditions
                colors = {}
                for i, condition in enumerate(legend_labels):
                    colors[condition] = color[i]
                    
                # Assign a color to each sample based on its condition
                sample_colors = [colors[condition] for sample in mds_results_df.index for condition, samples_list in conditions.items() if sample in samples_list]
                      
                # Scatter plot with samples colored by condition
                plt.figure(figsize=(12, 8))
                for condition, color in colors.items():
                    indices = [i for i, sample in enumerate(mds_results_df.index) if sample in conditions[condition]]
                    plt.scatter(mds_results_df.iloc[indices, 0], mds_results_df.iloc[indices, 1], color=color, label=condition)
                    # Add ellipses for each condition
                    if len(indices) > 1:  # Need at least 2 points to fit an ellipse
                        x_coords = mds_results_df.iloc[indices, 0]
                        y_coords = mds_results_df.iloc[indices, 1]
                        
                        # Calculate covariance matrix and mean
                        cov = np.cov(x_coords, y_coords)
                        print(cov)
                        eigenvalues, eigenvectors = np.linalg.eig(cov)
                        eigenvalues = np.sqrt(eigenvalues)
                        
                        # Calculate angle of ellipse
                        angle = np.rad2deg(np.arctan2(*eigenvectors[:, 0][::-1]))
                        
                        # Create an ellipse
                        ellipse = Ellipse(xy=(np.mean(x_coords), np.mean(y_coords)),
                                          width=eigenvalues[0] * 2, height=eigenvalues[1] * 2,
                                          angle=angle, edgecolor=color, facecolor=color, alpha=0.2, lw=2)
                        
                        plt.gca().add_patch(ellipse)
                    
                    
                plt.legend(title='Conditions', loc='best', bbox_to_anchor=(1, 1))
                plt.xlabel('NMDS 1')
                plt.ylabel('NMDS 2')
                plt.title('NMDS Plot')
                plt.annotate(f'Stress: {stress:.4f}', xy=(0.83, -0.06), xycoords='axes fraction', fontsize=12,
                 verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.5)) # Add an annotation for the stress value
                format_file = input("Which file format would you like to save the plot ? SVG or PDF or PNG")
                if format_file == 'PDF':
                    plt.savefig("NMDS_one_condition.pdf", format='pdf', bbox_inches='tight', pad_inches=0.2)
                elif format_file == 'SVG':
                    plt.savefig("NMDS_one_condition.svg", format='svg', bbox_inches='tight', pad_inches=0.2)
                elif format_file == 'PNG':
                    plt.savefig("NMDS_one_condition.png", format='png', bbox_inches='tight', pad_inches=0.2)
                plt.show()
            else :
                print("The stress variable is ", stress, ". It's greater than 0.2. Perform a PCoA analysis instead.")
                beta_diversity_df=beta_diversity.to_data_frame()
                # Replace row and column numbers with sample names
                sample_names = list(asv_info.columns[1:])
                beta_diversity_df.index = sample_names
                beta_diversity_df.columns = sample_names
                conditions=one_condition_struct(asv_info, condition) # Retrieve the structure containing sample data based on conditions
                
                # Perform PCoA on beta diversity
                pcoa_results = pcoa(beta_diversity)
                pcoa_results.samples.index = sample_names
                
                # Extract legend labels from conditions
                legend_labels = []
                for condition,samples_list in conditions.items():
                    legend_labels.append(condition)
                
                # Make unique colors from seaborn palette rainbow for each condition
                color = sns.color_palette('rainbow', n_colors=len(legend_labels)) # Generate as many colors as conditions
                colors = {}
                for i, condition in enumerate(legend_labels):
                    colors[condition] = color[i]
                    
                # Assign a color to each sample based on its condition
                sample_colors = [colors[condition] for sample in pcoa_results.samples.index for condition, samples_list in conditions.items() if sample in samples_list]
                      
                # Scatter plot with samples colored by condition
                plt.figure(figsize=(12, 8))
                for condition, color in colors.items():
                    indices = [i for i, sample in enumerate(pcoa_results.samples.index) if sample in conditions[condition]] # Find indices of samples belonging to the current condition
                    plt.scatter(pcoa_results.samples.iloc[indices, 0], pcoa_results.samples.iloc[indices, 1], color=color, label=condition)  # Scatter plot the samples belonging to the current condition using PC1 and PC2 coordinates
                    # Add ellipses for each condition
                    if len(indices) > 1:  # Need at least 2 points to fit an ellipse
                        x_coords = pcoa_results.samples.iloc[indices, 0]
                        y_coords = pcoa_results.samples.iloc[indices, 1]
                        
                        # Calculate covariance matrix and mean
                        cov = np.cov(x_coords, y_coords)
                        print(cov)
                        eigenvalues, eigenvectors = np.linalg.eig(cov)
                        eigenvalues = np.sqrt(eigenvalues)
                        
                        # Calculate angle of ellipse
                        angle = np.rad2deg(np.arctan2(*eigenvectors[:, 0][::-1]))
                        
                        # Create an ellipse
                        ellipse = Ellipse(xy=(np.mean(x_coords), np.mean(y_coords)),
                                          width=eigenvalues[0] * 2, height=eigenvalues[1] * 2,
                                          angle=angle, edgecolor=color, facecolor=color, alpha=0.2, lw=2)
                        
                        plt.gca().add_patch(ellipse)
                    
                plt.legend(title='Conditions', loc='best', bbox_to_anchor=(1, 1))
                plt.tight_layout()
                plt.title('PCoA Plot')
                plt.xlabel('Axis 1')
                plt.ylabel('Axis 2')
                format_file = input("Which file format would you like to save the plot ? SVG or PDF or PNG")
                if format_file == 'PDF':
                    plt.savefig("PCoA_one_condition.pdf", format='pdf', bbox_inches='tight', pad_inches=0.2)
                elif format_file == 'SVG':
                    plt.savefig("PCoA_one_condition.svg", format='svg', bbox_inches='tight', pad_inches=0.2)
                elif format_file == 'PNG':
                    plt.savefig("PCoA_one_condition.png", format='png', bbox_inches='tight', pad_inches=0.2)
                plt.show()
        elif number_condition == 'two':
            if stress <= 0.2:
                mds_results_df = pd.DataFrame(mds_results, columns=['Dimension 1', 'Dimension 2'])
                # Replace row and column numbers with sample names
                sample_names = list(asv_info.columns[1:])
                mds_results_df.index = sample_names
                    
                # Retrieve the structure containing sample data based on conditions
                conditions, cond1_name, cond2_name=two_conditions_struct(asv_info, condition) # Get the list of samples associated with the condition
                
                # Extract legend labels from conditions
                legend_labels1 = sorted(set(condition1 for condition1, _ in conditions.keys()))
                legend_labels2 = sorted(set(condition2 for _, condition2 in conditions.keys()))
                    
                # Make unique colors from seaborn palette rainbow for each value of the first condition
                color_palette = sns.color_palette('rainbow', n_colors=len(legend_labels1))
                colors = {label: color for label, color in zip(legend_labels1, color_palette)}
                print(colors)

                # Define markers for the second condition
                markers = ['o', 's', '^', 'D', 'v', '*', '<', '>', 'p', 'h', 'H', '8', 'd', '1', '2', '3', '4', '8', '+', 'x', 'X', '|', '_']  # Different marker styles
                marker_styles = {label: marker for label, marker in zip(legend_labels2, markers)}
                print(marker_styles)

                # Assign colors and markers to each sample based on its conditions
                sample_colors = [colors[condition1] for sample in mds_results_df.index for (condition1, condition2), samples_list in conditions.items() if sample in samples_list]
                sample_markers = [marker_styles[condition2] for sample in mds_results_df.index for (condition1, condition2), samples_list in conditions.items() if sample in samples_list]


                # Scatter plot with samples colored by the first condition and shaped by the second condition
                plt.figure(figsize=(12, 8))
                
                for cond1 in legend_labels1:
                    for cond2 in legend_labels2:
                        indices = [i for i, sample in enumerate(mds_results_df.index) if sample in conditions.get((cond1, cond2), [])]
                        if indices:
                            plt.scatter(mds_results_df.iloc[indices, 0], mds_results_df.iloc[indices, 1], color=[colors[cond1]]*len(indices), marker=marker_styles[cond2], label=f'{cond1}, {cond2}')
            
                # Create custom handles for the legend
                handles_colors = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=colors[cond], markersize=10, label=cond) for cond in legend_labels1]
                handles_markers = [plt.Line2D([0], [0], marker=marker_styles[cond], color='w', markerfacecolor='grey', markersize=10, label=cond) for cond in legend_labels2]
                
                # Add text labels for condition names
                handles = ([plt.Line2D([], [], color='none', label=cond1_name)] + handles_colors + [plt.Line2D([], [], color='none', label=cond2_name)] + handles_markers)
                labels = ([cond1_name] + legend_labels1 + [cond2_name] + legend_labels2)

                plt.legend(handles=handles, labels=labels, title='Conditions : ', loc='best', bbox_to_anchor=(1, 1))
                plt.tight_layout()
                plt.xlabel('NMDS 1')
                plt.ylabel('NMDS 2')
                plt.title('NMDS Plot')
                plt.annotate(f'Stress: {stress:.4f}', xy=(0.83, -0.06), xycoords='axes fraction', fontsize=12,
                 verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.5)) # Add an annotation for the stress value
                format_file = input("Which file format would you like to save the plot ? SVG or PDF or PNG")
                if format_file == 'PDF':
                    plt.savefig("NMDS_two_condition.pdf", format='pdf', bbox_inches='tight', pad_inches=0.2)
                elif format_file == 'SVG':
                    plt.savefig("NMDS_two_condition.svg", format='svg', bbox_inches='tight', pad_inches=0.2)
                elif format_file == 'PNG':
                    plt.savefig("NMDS_two_condition.png", format='png', bbox_inches='tight', pad_inches=0.2)
                plt.show()
            else :
                print("The stress variable is ", stress, ". It's greater than 0.2. Perform a PCoA analysis instead.")
                beta_diversity_df=beta_diversity.to_data_frame()
                # Replace row and column numbers with sample names
                sample_names = list(asv_info.columns[1:])
                beta_diversity_df.index = sample_names
                beta_diversity_df.columns = sample_names
                conditions, cond1_name, cond2_name =two_conditions_struct(asv_info, condition) # Get the list of samples associated with the condition
                
                # Perform PCoA on beta diversity
                pcoa_results = pcoa(beta_diversity)
                pcoa_results.samples.index = sample_names

                # Extract legend labels from conditions
                legend_labels1 = sorted(set(condition1 for condition1, _ in conditions.keys()))
                legend_labels2 = sorted(set(condition2 for _, condition2 in conditions.keys()))

                # Make unique colors from seaborn palette rainbow for each value of the first condition
                color_palette = sns.color_palette('rainbow', n_colors=len(legend_labels1))
                colors = {label: color for label, color in zip(legend_labels1, color_palette)}
                print(colors)

                # Define markers for the second condition
                markers = ['o', 's', '^', 'D', 'v', '*', '<', '>', 'p', 'h', 'H', '8', 'd', '1', '2', '3', '4', '8', '+', 'x', 'X', '|', '_']  # Different marker styles
                marker_styles = {label: marker for label, marker in zip(legend_labels2, markers)}
                print(marker_styles)

                # Assign colors and markers to each sample based on its conditions
                sample_colors = [colors[condition1] for sample in pcoa_results.samples.index for (condition1, condition2), samples_list in conditions.items() if sample in samples_list]
                sample_markers = [marker_styles[condition2] for sample in pcoa_results.samples.index for (condition1, condition2), samples_list in conditions.items() if sample in samples_list]

                # Scatter plot with samples colored by the first condition and shaped by the second condition
                plt.figure(figsize=(12, 8))
                    
                for cond1 in legend_labels1:
                    for cond2 in legend_labels2:
                        indices = [i for i, sample in enumerate(pcoa_results.samples.index) if sample in conditions.get((cond1, cond2), [])]
                        if indices:
                            plt.scatter(pcoa_results.samples.iloc[indices, 0], pcoa_results.samples.iloc[indices, 1], color=[colors[cond1]]*len(indices), marker=marker_styles[cond2], label=f'{cond1}, {cond2}')
                
                # Create custom handles for the legend
                handles_colors = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=colors[cond], markersize=10, label=cond) for cond in legend_labels1]
                handles_markers = [plt.Line2D([0], [0], marker=marker_styles[cond], color='w', markerfacecolor='grey', markersize=10, label=cond) for cond in legend_labels2]
                
                # Add text labels for condition names
                handles = ([plt.Line2D([], [], color='none', label=cond1_name)] + handles_colors + [plt.Line2D([], [], color='none', label=cond2_name)] + handles_markers)
                labels = ([cond1_name] + legend_labels1 + [cond2_name] + legend_labels2)

                plt.legend(handles=handles, labels=labels, title='Conditions : ', loc='best', bbox_to_anchor=(1, 1))
                plt.tight_layout()
                plt.title('PCoA Plot')
                plt.xlabel('Axis 1')
                plt.ylabel('Axis 2')
                format_file = input("Which file format would you like to save the plot? SVG or PDF or PNG")
                if format_file == 'PDF':
                    plt.savefig("PCoA_two_condition.pdf", format='pdf', bbox_inches='tight', pad_inches=0.2)
                elif format_file == 'SVG':
                    plt.savefig("PCoA_two_condition.svg", format='svg', bbox_inches='tight', pad_inches=0.2)
                elif format_file == 'PNG':
                    plt.savefig("PCoA_two_condition.png", format='png', bbox_inches='tight', pad_inches=0.2)
                plt.show()
        
        
    elif beta_representation == 'PCoA': # Create a PCoA plot of beta diversity, colored by samples
        pcoa_results = pcoa(beta_diversity) # Perform PCoA
        sample_names = list(asv_info.columns[1:]) # Extract sample name
        pcoa_results.samples.index = sample_names # Set the IDs of PCoA results samples to the sample names
        colors = sns.color_palette('rainbow', n_colors=len(sample_names)) # Generate as many colors as samples
        legend_handles = []
        plt.figure(figsize=(12, 8))
        
        for sample, color in zip(sample_names, colors): # Iterate over sample names and corresponding colors
           x = pcoa_results.samples.loc[sample, 'PC1']
           y = pcoa_results.samples.loc[sample, 'PC2']
           handle = plt.scatter(x, y, color=color, label=sample)
           legend_handles.append(handle)
            
        plt.legend(handles=legend_handles, title='Samples', loc='best', bbox_to_anchor=(1, 1))
        plt.tight_layout()
        plt.title('PCoA Plot')
        plt.xlabel('Axis 1')
        plt.ylabel('Axis 2')
        format_file = input("Which file format would you like to save the plot ? SVG or PDF or PNG")
        if format_file == 'PDF':
            plt.savefig("PCoA.pdf", format='pdf', bbox_inches='tight', pad_inches=0.2)
        elif format_file == 'SVG':
            plt.savefig("PCoA.svg", format='svg', bbox_inches='tight', pad_inches=0.2)
        elif format_file == 'PNG':
            plt.savefig("PCoA.png", format='png', bbox_inches='tight', pad_inches=0.2)
        plt.show()
    
    elif beta_representation == 'PCoA_grouped_by_condition': # Create a PCoA plot of beta diversity grouped by conditions
        beta_diversity_df=beta_diversity.to_data_frame()
        # Replace row and column numbers with sample names
        sample_names = list(asv_info.columns[1:])
        beta_diversity_df.index = sample_names
        beta_diversity_df.columns = sample_names
        
        number_condition = input("How many conditions do you have? (one / two): ")
        if number_condition == 'one':
            conditions=one_condition_struct(asv_info, condition) # Retrieve the structure containing sample data based on conditions
            
            # Perform PCoA on beta diversity
            pcoa_results = pcoa(beta_diversity)
            pcoa_results.samples.index = sample_names
            
            # Extract legend labels from conditions
            legend_labels = []
            for condition,samples_list in conditions.items():
                legend_labels.append(condition)
            
            # Make unique colors from seaborn palette rainbow for each condition
            color = sns.color_palette('rainbow', n_colors=len(legend_labels)) # Generate as many colors as conditions
            colors = {}
            for i, condition in enumerate(legend_labels):
                colors[condition] = color[i]
                
            # Assign a color to each sample based on its condition
            sample_colors = [colors[condition] for sample in pcoa_results.samples.index for condition, samples_list in conditions.items() if sample in samples_list]
            
            # Determine which axes to plot
            x_axis = int(input("Which PCoA axis do you want to display on x-axis ( axis 1 = 0 / axis 2 = 1 / axis 3 = 2) ?"))
            y_axis = int(input("Which PCoA axis do you want to display on y-axis ( axis 1 = 0 / axis 2 = 1 / axis 3 = 2) ?"))

            # Scatter plot with samples colored by condition
            plt.figure(figsize=(12, 8))
            for condition, color in colors.items():
                indices = [i for i, sample in enumerate(pcoa_results.samples.index) if sample in conditions[condition]] # Find indices of samples belonging to the current condition
                plt.scatter(pcoa_results.samples.iloc[indices, x_axis], pcoa_results.samples.iloc[indices, y_axis], color=color, label=condition)  # Scatter plot the samples belonging to the current condition using PC1 and PC2 coordinates
                # Add ellipses for each condition
                if len(indices) > 1:  # Need at least 2 points to fit an ellipse
                    x_coords = pcoa_results.samples.iloc[indices, x_axis]
                    y_coords = pcoa_results.samples.iloc[indices, y_axis]
                    
                    # Calculate covariance matrix and mean
                    cov = np.cov(x_coords, y_coords)
                    eigenvalues, eigenvectors = np.linalg.eig(cov)
                    eigenvalues = np.sqrt(eigenvalues)
                    
                    # Calculate angle of ellipse
                    angle = np.rad2deg(np.arctan2(*eigenvectors[:, 0][::-1]))
                    
                    # Create an ellipse
                    ellipse = Ellipse(xy=(np.mean(x_coords), np.mean(y_coords)),
                                      width=eigenvalues[0] * 2, height=eigenvalues[1] * 2,
                                      angle=angle, edgecolor=color, facecolor=color, alpha=0.2, lw=2)
                    
                    plt.gca().add_patch(ellipse)
            
            plt.legend(title='Conditions', loc='best', bbox_to_anchor=(1, 1))
            plt.tight_layout()
            plt.title('PCoA Plot')
            plt.xlabel(f'Axis {x_axis + 1}')
            plt.ylabel(f'Axis {y_axis + 1}')
            format_file = input("Which file format would you like to save the plot ? SVG or PDF or PNG")
            if format_file == 'PDF':
                plt.savefig("PCoA_one_condition.pdf", format='pdf', bbox_inches='tight', pad_inches=0.2)
            elif format_file == 'SVG':
                plt.savefig("PCoA_one_condition.svg", format='svg', bbox_inches='tight', pad_inches=0.2)
            elif format_file == 'PNG':
                plt.savefig("PCoA_one_condition.png", format='png', bbox_inches='tight', pad_inches=0.2)
            plt.show()
            
        elif number_condition == 'two':
            conditions, cond1_name, cond2_name =two_conditions_struct(asv_info, condition) # Get the list of samples associated with the condition
            
            # Perform PCoA on beta diversity
            pcoa_results = pcoa(beta_diversity)
            pcoa_results.samples.index = sample_names

            # Extract legend labels from conditions
            legend_labels1 = sorted(set(condition1 for condition1, _ in conditions.keys()))
            legend_labels2 = sorted(set(condition2 for _, condition2 in conditions.keys()))

            # Make unique colors from seaborn palette rainbow for each value of the first condition
            color_palette = sns.color_palette('rainbow', n_colors=len(legend_labels1))
            colors = {label: color for label, color in zip(legend_labels1, color_palette)}
            print(colors)

            # Define markers for the second condition
            markers = ['o', 's', '^', 'D', 'v', '*', '<', '>', 'p', 'h', 'H', '8', 'd', '1', '2', '3', '4', '8', '+', 'x', 'X', '|', '_']  # Different marker styles
            marker_styles = {label: marker for label, marker in zip(legend_labels2, markers)}
            print(marker_styles)

            # Assign colors and markers to each sample based on its conditions
            sample_colors = [colors[condition1] for sample in pcoa_results.samples.index for (condition1, condition2), samples_list in conditions.items() if sample in samples_list]
            sample_markers = [marker_styles[condition2] for sample in pcoa_results.samples.index for (condition1, condition2), samples_list in conditions.items() if sample in samples_list]

            # Scatter plot with samples colored by the first condition and shaped by the second condition
            plt.figure(figsize=(12, 8))
                
            for cond1 in legend_labels1:
                for cond2 in legend_labels2:
                    indices = [i for i, sample in enumerate(pcoa_results.samples.index) if sample in conditions.get((cond1, cond2), [])]
                    if indices:
                        plt.scatter(pcoa_results.samples.iloc[indices, 0], pcoa_results.samples.iloc[indices, 1], color=[colors[cond1]]*len(indices), marker=marker_styles[cond2], label=f'{cond1}, {cond2}')
            
            # Create custom handles for the legend
            handles_colors = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=colors[cond], markersize=10, label=cond) for cond in legend_labels1]
            handles_markers = [plt.Line2D([0], [0], marker=marker_styles[cond], color='w', markerfacecolor='grey', markersize=10, label=cond) for cond in legend_labels2]
            
            # Add text labels for condition names
            handles = ([plt.Line2D([], [], color='none', label=cond1_name)] + handles_colors + [plt.Line2D([], [], color='none', label=cond2_name)] + handles_markers)
            labels = ([cond1_name] + legend_labels1 + [cond2_name] + legend_labels2)

            plt.legend(handles=handles, labels=labels, title='Conditions : ', loc='best', bbox_to_anchor=(1, 1))
            plt.tight_layout()
            plt.title('PCoA Plot')
            plt.xlabel('Axis 1')
            plt.ylabel('Axis 2')
            format_file = input("Which file format would you like to save the plot? SVG or PDF or PNG")
            if format_file == 'PDF':
                plt.savefig("PCoA_two_condition.pdf", format='pdf', bbox_inches='tight', pad_inches=0.2)
            elif format_file == 'SVG':
                plt.savefig("PCoA_two_condition.svg", format='svg', bbox_inches='tight', pad_inches=0.2)
            elif format_file == 'PNG':
                plt.savefig("PCoA_two_condition.png", format='png', bbox_inches='tight', pad_inches=0.2)
            plt.show()
        

    else:
        print("Beta diversity representation not supported.")
        exit()
    


#################################
## Graph taxonomic composition ##
#################################


def get_taxo_count(database, taxo_level, samples, number_sample): # Function to handles missing values by replacing them with information from a higher taxonomic level
    if number_sample == 1:
        samples = [samples]
    else:
        samples = list(samples)

    # Fonction to replace NaN values with the higher level taxon
    def replace_na(row, taxo_level_index):
        for i, level in enumerate(taxonomy_levels[taxo_level_index-1::-1]): # Iterate backwards through the taxonomic levels until a non-NaN value is found
            if level in row and pd.notna(row[level]):
                return f"{row[level]} ({taxonomy_levels[taxo_level_index - 1 - i]})" # Return the taxon with the level in parentheses
        return row[taxo_level]

    taxo_level_index = taxonomy_levels.index(taxo_level)  # Get the index of the selected taxonomic level
    relevant_taxonomy_levels = taxonomy_levels[:taxo_level_index + 1] # Create a list of relevant taxonomic levels up to the selected level
    taxo_count = database[relevant_taxonomy_levels + samples] # Extract the relevant taxonomic levels and sample counts from the file database

    # Apply the replace_na function to each row to replace NaN values with the higher level taxon
    taxo_count[taxo_level] = taxo_count.apply(
        lambda row: replace_na(row, taxo_level_index) if pd.isna(row[taxo_level]) else row[taxo_level], axis=1
        )

    taxo_count = taxo_count[[taxo_level] + samples] # Remove columns corresponding to higher taxonomic levels
    return taxo_count

def barplot_one_sample(database): # Function to create a bar plot of taxonomic composition for a given sample
    print("Available samples :")
    available_samples = list(database.columns[database.columns.get_loc('sequence')+1:])
    for sample in available_samples:
        print(sample)
    sample = input("Which sample would you like to analyze ? ")
    taxo_level = input(f'What taxonomic level do you want ? \n {taxonomy_levels} \nYour choice : ')
    
    taxo_count = get_taxo_count(database, taxo_level, sample, 1)
    print(taxo_count.tail(50))
    taxo_sum = taxo_count.groupby(taxo_level).sum() # Group by taxonomic level and sum counts
    taxo_sum[sample] = taxo_sum[sample].astype(float) # Convert count values to float
    taxo_proportions = taxo_sum / taxo_sum.sum() # Calculate proportions of each taxonomic category
    
    # Selection of the x top most represented categories :
    nb_taxa = int(input("How many taxa do you want to display ? "))
    top_x_taxa = taxo_proportions.apply(lambda x: x.nlargest(nb_taxa), axis=0) # Keep the first x categories
    other_taxa = 1 - top_x_taxa.sum(axis=0) # Calculation of the sum of the proportions of the remaining categories for each sample
    print(other_taxa)
    if (other_taxa > 0).any():
        top_x_taxa.loc['Others'] = other_taxa.values # Addition of the 'Others' line with the proportions of the 'Others' category for each sample

    # Bar plot of taxonomic proportions
    plt.figure(figsize=(10, 6))
    color = sns.color_palette('rainbow', n_colors=nb_taxa+1) # Generate as many colors as categories
    top_x_taxa[sample].sort_values(ascending=False).plot(color=color, kind='bar')
    plt.xlabel(taxo_level)
    plt.ylabel('Proportion')
    plt.title(f'Proportion of {taxo_level} taxonomic category in sample {sample}')
    plt.xticks(rotation=90, ha='right')
    plt.tight_layout()
    format_file = input("Which file format would you like to save the plot ? SVG or PDF or PNG")
    if format_file == 'PDF':
        plt.savefig("barplot_one_sample.pdf", format='pdf', pad_inches=0.2)
    elif format_file == 'SVG':
        plt.savefig("barplot_one_sample.svg", format='svg', pad_inches=0.2)
    elif format_file == 'PNG':
        plt.savefig("barplot_one_sample.png", format='png', pad_inches=0.2)
    plt.show()

    
def barplot_all_samples(database): # Function to create a bar plot of taxonomic composition for all samples
    # Select columns corresponding to samples (after the "sequence" column)
    sequence_index = list(database.columns).index("sequence")
    samples_columns = database.columns[sequence_index + 1:]
    taxo_level = input(f'What taxonomic level do you want ? \n {taxonomy_levels} \nYour choice : ')
    
    taxo_count = get_taxo_count(database, taxo_level, samples_columns, 2) # Extract taxonomic counts for the selected taxonomic level and all samples
    taxo_sum = taxo_count.groupby(taxo_level).sum() # Group by taxonomic level and sum counts
    taxo_proportions = taxo_sum.div(taxo_sum.sum(axis=0), axis=1) # Calculate proportions of each taxonomic category for each sample

    # Selection of the x most represented categories for each sample :
    nb_taxa = int(input("How many taxa do you want to display ? "))
    total_abundance = taxo_proportions.sum(axis=1) # Calculate the total abundance of each taxon across all samples
    top_x_abundance = total_abundance.nlargest(nb_taxa).index  # Get the top x most abundant taxons
    top_x_taxa = taxo_proportions.loc[top_x_abundance] # Keep only the rows corresponding to the top x taxa
    other_taxa = 1 - top_x_taxa.sum(axis=0) # Calculation of the sum of the proportions of the remaining categories for each sample
    if (other_taxa > 0).any():
        top_x_taxa.loc['Others'] = other_taxa.values # Addition of the 'Others' line with the proportions of the 'Others' category for each sample
    top_x_taxa_transpose = top_x_taxa.transpose() # Transposition to have samples in rows and categories in columns
        
    # Bar plot of taxonomic proportions for all samples
    fig, ax = plt.subplots(figsize=(18, 12))
    color = sns.color_palette('rainbow', n_colors=nb_taxa+1) # Generate as many colors as categories
    top_x_taxa_transpose.plot(kind='bar', stacked=True, color=color, ax=ax)
    plt.xlabel('Samples')
    plt.ylabel('Proportion')
    plt.title(f'Proportion of taxonomic category {taxo_level} for all samples')
    plt.legend(title=taxo_level, loc='best', bbox_to_anchor=(1, 1))
    plt.xticks(rotation=90, ha='right')
    plt.tight_layout()
    format_file = input("Which file format would you like to save the plot ? SVG or PDF or PNG")
    if format_file == 'PDF':
        plt.savefig("barplot_composition_all_samples.pdf", format='pdf', bbox_inches='tight', pad_inches=0.2)
    elif format_file == 'SVG':
        plt.savefig("barplot_composition_all_samples.svg", format='svg', bbox_inches='tight', pad_inches=0.2)
    elif format_file == 'PNG':
        plt.savefig("barplot_composition_all_samples.png", format='png', bbox_inches='tight', pad_inches=0.2)
    plt.show()
        
        

def barplot_all_samples_condition(database, condition, asv_info):
    # Extract column names of samples
    sequence_index = list(database.columns).index("sequence")
    samples_columns = database.columns[sequence_index + 1:]
    
    taxo_level = input(f'What taxonomic level do you want ? \n {taxonomy_levels} \nYour choice : ')
    
    conditions=one_condition_struct(asv_info, condition) # Get the list of samples associated with the condition
    
    plot_level = input("Do you want a merged plot with all conditions (merged_plot) or separate plots for each condition (separate_plots) or plot with average proportion for each categories (total_abundance) ? ")
    
    if plot_level == 'total_abundance':
        # Initialize a dictionary to store summed abundances for each condition
        condition_abundances = {cond: 0 for cond in conditions}
        
        # Iterate over each condition
        for i, (condition, samples_list) in enumerate(conditions.items()):
            # Select samples corresponding to the condition
            taxo_count = get_taxo_count(database, taxo_level, samples_list, 2) # Extract taxonomic counts for the selected taxonomic level and all samples
            taxo_sum = taxo_count.groupby(taxo_level).sum()
            taxo_proportions = taxo_sum.div(taxo_sum.sum(axis=0), axis=1)
            taxo_proportions_transpose = taxo_proportions.transpose()
            condition_abundances[condition] = taxo_proportions.sum(axis=1)/len(samples_list) # Sum proportions across samples for each taxonomic category
                
        # Plot the barplot
        plt.figure(figsize=(12, 8))
        legend_handles = []  # List to store legend handles
        
        all_taxons = set()
        for taxon_abundances in condition_abundances.values():
            all_taxons.update(taxon_abundances.index)

        #color_map = sns.color_palette("rainbow", len(all_taxons))
        
        nb_taxa = int(input("How many taxa do you want to display ? "))
        color_map = sns.color_palette('rainbow', n_colors=nb_taxa+1) # Generate as many colors as categories
        for i, (condition_name, taxon_abundances) in enumerate(condition_abundances.items()):
            top_x_taxa = taxon_abundances.sort_values(ascending=False)[:nb_taxa] # Keep the first x categories
            other_taxa = 1 - top_x_taxa.sum(axis=0)
            if (other_taxa > 0).any():
                top_x_taxa.loc['Others'] = other_taxa # Addition of the 'Others' line with the proportions of the 'Others' category for each condition
            
            # Initialize y_offset to keep track of the vertical position of the next bar
            y_offset = 0
            for j, (taxon, abundance) in enumerate(top_x_taxa.items()):
                color = color_map[j] # Assign color based on taxon
                plt.bar(condition_name, abundance, color=color, bottom=y_offset)  # Plot each taxon with assigned color
                y_offset += abundance
                if i == 0:  # Only add handles from the first condition to avoid duplicates
                    legend_handles.append(taxon)
            
        plt.xlabel('Condition')
        plt.ylabel('Total Proportion')
        plt.title('Total Proportion by Condition')
        plt.xticks(rotation=90)
        plt.legend(legend_handles, title=taxo_level, loc='best', bbox_to_anchor=(1, 1))
        plt.tight_layout()
        
        format_file = input("Which file format would you like to save the plot ? SVG or PDF or PNG")
        if format_file == 'PDF':
            plt.savefig("barplot_total_abundance_by_condition.pdf", format='pdf', bbox_inches='tight', pad_inches=0.2)
        elif format_file == 'SVG':
            plt.savefig("barplot_total_abundance_by_condition.svg", format='svg', bbox_inches='tight', pad_inches=0.2)
        elif format_file == 'PNG':
            plt.savefig("barplot_total_abundance_by_condition.png", format='png', bbox_inches='tight', pad_inches=0.2)
        plt.show()
        
    if plot_level == 'merged_plot':
        fig, axs = plt.subplots(1, len(conditions), figsize=(18, 12), sharey=True)
        num_conditions = len(conditions)
        num_cols = len(conditions) #2  # Nombre de colonnes dans la grille
        num_rows = (num_conditions + num_cols - 1) // num_cols  # Calcul du nombre de lignes nécessaire
        
        taxo_count = get_taxo_count(database, taxo_level, samples_columns, 2) # Extract taxonomic counts for the selected taxonomic level and all samples
        taxo_sum = taxo_count.groupby(taxo_level).sum() # Group by taxonomic level and sum counts
        taxo_proportions = taxo_sum.div(taxo_sum.sum(axis=0), axis=1) # Calculate proportions of each taxonomic category for each sample
        
        # Select top x taxonomic categories
        nb_taxa = int(input("How many taxa do you want to display ? "))
        color = sns.color_palette('rainbow', n_colors=nb_taxa+1) # Generate as many colors as categories
        total_abundance = taxo_proportions.sum(axis=1) # Calculate the total abundance of each taxon across all samples
        top_x_abundance = total_abundance.nlargest(nb_taxa).index  # Get the top x most abundant taxons
        top_x_taxa = taxo_proportions.loc[top_x_abundance] # Keep only the rows corresponding to the top x taxa
        other_taxa = 1 - top_x_taxa.sum(axis=0) # Calculation of the sum of the proportions of the remaining categories for each sample
        if (other_taxa > 0).any():
            top_x_taxa.loc['Others'] = other_taxa.values # Addition of the 'Others' line with the proportions of the 'Others' category for each sample
        top_x_taxa_transpose = top_x_taxa.transpose() # Transposition to have samples in rows and categories in columns
        
        # Associate one color to each taxa
        taxon_colors = {} # Dictionary to store unique colors for each taxon
        for i, (taxon, sample) in enumerate(top_x_taxa_transpose.items()):
            taxon_colors[taxon] = color[i]  # Associer le taxon à la couleur générée
        
        # Barplots for each condition
        for i, (condition, samples_list) in enumerate(conditions.items()):
            row = i // num_cols  # Ligne actuelle dans la grille
            col = i % num_cols   # Colonne actuelle dans la grille
        
            print(condition)
            print(samples_list)
        
            # Plot the barplot for the current condition
            ax = axs[i]
            top_x_taxa_transpose.loc[samples_list].plot(kind='bar', stacked=True, color=[taxon_colors[taxon] for taxon in top_x_taxa_transpose.columns], ax=ax, label=condition, legend=None)
            ax.set_xlabel('Samples')
            ax.set_ylabel('Abundance')
            ax.set_title(condition)
            ax.tick_params(axis='x', rotation=90)
        
        plt.tight_layout()
        # Créer une légende commune
        legend_elements = [plt.Line2D([0], [0], color=taxon_colors[taxon], lw=4, label=taxon) for taxon in taxon_colors]
        plt.legend(handles=legend_elements, title=taxo_level, loc='best', bbox_to_anchor=(1, 1))

        # Sauvegarder et afficher le plot
        format_file = input("Which file format would you like to save the plot ? SVG or PDF or PNG")
        if format_file == 'PDF':
            plt.savefig("barplot_merged_conditions.pdf", format='pdf', bbox_inches='tight', pad_inches=0.2)
        elif format_file == 'SVG':
            plt.savefig("barplot_merged_conditions.svg", format='svg', bbox_inches='tight', pad_inches=0.2)
        elif format_file == 'PNG':
            plt.savefig("barplot_merged_conditions.png", format='png', bbox_inches='tight', pad_inches=0.2)
        plt.show()


    elif plot_level == 'separate_plots':
        # Barplots for each condition
        for condition, samples_list in conditions.items():
            taxo_count = get_taxo_count(database, taxo_level, samples_list, 2) # Extract taxonomic counts for the selected taxonomic level and all samples
            taxo_sum = taxo_count.groupby(taxo_level).sum()
            taxo_proportions = taxo_sum.div(taxo_sum.sum(axis=0), axis=1)

            # Select top x taxonomic categories
            nb_taxa = int(input("How many taxa do you want to display ? "))
            color = sns.color_palette('rainbow', n_colors=nb_taxa+1) # Generate as many colors as categories
            top_x_taxa = taxo_proportions.sum(axis=1).nlargest(nb_taxa)
            taxo_proportions = taxo_proportions.transpose() # Transposition to have samples in rows and categories in columns
            top_x_taxa = taxo_proportions[top_x_taxa.index]
            other_taxa = 1 - top_x_taxa.sum(axis=1) # Calculation of the sum of the proportions of the remaining categories for each sample
            top_x_taxa['Others'] = other_taxa.values # Addition of the 'Others' line with the proportions of the 'Others' category for each sample

            # Plot the barplot for the current condition
            fig, ax = plt.subplots(figsize=(18, 12))
            top_x_taxa.plot(kind='bar', stacked=True, color=color, ax=ax, label=condition)
            plt.xlabel('Samples')
            plt.ylabel('Abundance')
            plt.title(f'Abundance of taxonomic category {taxo_level} for condition : {condition}')
            plt.legend(title=taxo_level, loc='best', bbox_to_anchor=(1, 1))
            plt.xticks(rotation=90, ha='right')
            plt.tight_layout()
            format_file = input("Which file format would you like to save the plot ? SVG or PDF or PNG")
            if format_file == 'PDF':
                plt.savefig(f"barplot_{condition}.pdf", format='pdf', bbox_inches='tight', pad_inches=0.2)
            elif format_file == 'SVG':
                plt.savefig(f"barplot_{condition}.svg", format='svg', bbox_inches='tight', pad_inches=0.2)
            elif format_file == 'PNG':
                plt.savefig(f"barplot_{condition}.png", format='png', bbox_inches='tight', pad_inches=0.2)
            plt.show()
    else:
        print("Plots representation not supported.")
        
        

def barplot_all_samples_two_conditions(database, condition, asv_info):
    # Extract column names of samples
    sequence_index = list(database.columns).index("sequence")
    samples_columns = database.columns[sequence_index + 1:]
    
    taxo_level = input(f'What taxonomic level do you want ? \n {taxonomy_levels} \nYour choice : ')
    
    conditions=two_conditions_struct(asv_info, condition) # Get the list of samples associated with the condition
    
    plot_level = input("Do you want a plot with average proportion (total_abudance) or separate plots for each condition (separate_plots) ? ")
    
    if plot_level == 'total_abudance':
        # Initialize a dictionary to store summed abundances for each condition
        condition_abundances = {cond: 0 for cond in conditions}
        
        # Iterate over each conditions
        for i, ((cond1, cond2), samples_list) in enumerate(conditions.items()):
            # Select samples corresponding to the condition
            taxo_count = get_taxo_count(database, taxo_level, samples_list, 2) # Extract taxonomic counts for the selected taxonomic level and all samples
            taxo_sum = taxo_count.groupby(taxo_level).sum()
            taxo_proportions = taxo_sum.div(taxo_sum.sum(axis=0), axis=1)
            taxo_proportions_transpose = taxo_proportions.transpose()
            condition_abundances[(cond1, cond2)] = taxo_proportions.sum(axis=1) / len(samples_list)  # Somme des proportions pour chaque catégorie taxonomique
       
        # Plot the barplot
        plt.figure(figsize=(12, 8))
        legend_handles = []  # List to store legend handles
        
        all_taxons = set()
        for taxon_abundances in condition_abundances.values():
            all_taxons.update(taxon_abundances.index)
        
        nb_taxa = int(input("How many taxa do you want to display ? "))
        color_map = sns.color_palette('rainbow', n_colors=nb_taxa+1) # Generate as many colors as categories
        for i, (condition_name, taxon_abundances) in enumerate(condition_abundances.items()):
            top_x_taxa = taxon_abundances.sort_values(ascending=False)[:nb_taxa]
            other_taxa = 1 - top_x_taxa.sum(axis=0)
            if (other_taxa > 0).any():
                top_x_taxa.loc['Others'] = other_taxa # Addition of the 'Others' line with the proportions of the 'Others' category for each condition
            
            # Initialize y_offset to keep track of the vertical position of the next bar
            y_offset = 0
            for j, (taxon, abundance) in enumerate(top_x_taxa.items()):
                color = color_map[j] # Assign color based on taxon
                plt.bar((f"{condition_name[0]} - {condition_name[1]}"), abundance, color=color, bottom=y_offset)  # Plot each taxon with assigned color
                y_offset += abundance
                if i == 0:  # Only add handles from the first condition to avoid duplicates
                    legend_handles.append(taxon)
            
            plt.text(i, -0.05, condition_name[0], ha='center', va='bottom', fontsize=12, fontweight='bold')
            plt.text(i, 1.04, condition_name[1], ha='center', va='top', fontsize=12, fontweight='bold')

                    
        plt.xlabel('Condition')
        plt.ylabel('Total Proportion')
        plt.title('Total Proportion by Conditions')
        plt.xticks(rotation=90)
        plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        plt.legend(legend_handles, title=taxo_level, loc='best', bbox_to_anchor=(1, 1))
        plt.tight_layout()
        
        format_file = input("Which file format would you like to save the plot ? SVG or PDF or PNG")
        if format_file == 'PDF':
            plt.savefig("barplot_total_abundance_by_two_condition.pdf", format='pdf', bbox_inches='tight', pad_inches=0.2)
        elif format_file == 'SVG':
            plt.savefig("barplot_total_abundance_by_two_condition.svg", format='svg', bbox_inches='tight', pad_inches=0.2)
        elif format_file == 'PNG':
            plt.savefig("barplot_total_abundance_by_two_condition.png", format='png', bbox_inches='tight', pad_inches=0.2)
        plt.show()

    elif plot_level == 'separate_plots':
        # Barplots for each condition
        for condition, samples_list in conditions.items():
        
            taxo_count = get_taxo_count(database, taxo_level, samples_list, 2) # Extract taxonomic counts for the selected taxonomic level and all samples
            taxo_sum = taxo_count.groupby(taxo_level).sum()
            taxo_proportions = taxo_sum.div(taxo_sum.sum(axis=0), axis=1)

            # Select top x taxonomic categories
            nb_taxa = int(input("How many taxa do you want to display ? "))
            color = sns.color_palette('rainbow', n_colors=nb_taxa+1) # Generate as many colors as categories
            top_x_taxa = taxo_proportions.sum(axis=1).nlargest(nb_taxa)
            taxo_proportions = taxo_proportions.transpose() # Transposition to have samples in rows and categories in columns
            top_x_taxa = taxo_proportions[top_x_taxa.index]
            other_taxa = 1 - top_x_taxa.sum(axis=1) # Calculation of the sum of the proportions of the remaining categories for each sample
            top_x_taxa['Others'] = other_taxa.values # Addition of the 'Others' line with the proportions of the 'Others' category for each sample

            # Plot the barplot for the current condition
            fig, ax = plt.subplots(figsize=(18, 12))
            top_x_taxa.plot(kind='bar', stacked=True, color=color, ax=ax, label=condition)
            plt.xlabel('Samples')
            plt.ylabel('Abundance')
            plt.title(f'Abundance of taxonomic category {taxo_level} for condition: {condition}')
            plt.legend(title=taxo_level, loc='best', bbox_to_anchor=(1, 1))
            plt.xticks(rotation=90, ha='right')
            plt.tight_layout()
            format_file = input("Which file format would you like to save the plot ? SVG or PDF or PNG")
            if format_file == 'PDF':
                plt.savefig(f"barplot_composition_{condition}.pdf", format='pdf', bbox_inches='tight', pad_inches=0.2)
            elif format_file == 'SVG':
                plt.savefig(f"barplot_composition_{condition}.svg", format='svg', bbox_inches='tight', pad_inches=0.2)
            elif format_file == 'PNG':
                plt.savefig(f"barplot_composition_{condition}.png", format='png', bbox_inches='tight', pad_inches=0.2)
            plt.show()
    else:
        print("Plots representation not supported.")
        
        

def heatmap(database): # Function to create a heatmap of taxonomic composition for all samples
    # Extract column names of samples
    sequence_index = list(database.columns).index("sequence")
    samples_columns = database.columns[sequence_index + 1:]
    
    taxo_level = input(f'What taxonomic level do you want ? \n {taxonomy_levels} \nYour choice : ')
    
    # Calculate proportions of each taxonomic category for each sample
    taxo_count = get_taxo_count(database, taxo_level, samples_columns, 2) # Extract taxonomic counts for the selected taxonomic level and all samples
    taxo_sum = taxo_count.groupby(taxo_level).sum()
    taxo_proportions = taxo_sum.div(taxo_sum.sum(axis=0), axis=1)
    taxo_proportions_heatmap = taxo_proportions.transpose()

    # Heatmap of taxonomic proportions for all samples
    plt.figure(figsize=(12, 8))
    sns.heatmap(taxo_proportions_heatmap, cmap='rainbow')
    plt.xlabel(taxo_level)
    plt.ylabel('Samples')
    plt.title(f'Heatmap of the proportion of taxonomic category {taxo_level} in each sample')
    plt.xticks(rotation=90, ha='right')
    plt.tight_layout()
    format_file = input("Which file format would you like to save the plot ? SVG or PDF or PNG")
    if format_file == 'PDF':
        plt.savefig("heatmap_proportion.pdf", format='pdf', bbox_inches='tight', pad_inches=0.2)
    elif format_file == 'SVG':
        plt.savefig("heatmap_proportion.svg", format='svg', bbox_inches='tight', pad_inches=0.2)
    elif format_file == 'PNG':
        plt.savefig("heatmap_proportion.png", format='png', bbox_inches='tight', pad_inches=0.2)
    plt.show()


def piechart(database): # Function to create a pie chart of taxonomic composition for all samples
    # Extract column names of samples
    sequence_index = list(database.columns).index("sequence")
    samples_columns = database.columns[sequence_index + 1:]
    
    taxo_level = input(f'What taxonomic level do you want ? \n {taxonomy_levels} \nYour choice : ')
    
    nb_taxa = int(input("How many taxa do you want to display ? "))
    
    taxo_count = get_taxo_count(database, taxo_level, samples_columns, 2) # Extract taxonomic counts for the selected taxonomic level and all samples
    taxo_sum = taxo_count.groupby(taxo_level).sum() # Group by the chosen taxonomic level and sum the counts
    total_counts = taxo_sum.sum(axis=1) # Sum the counts across all samples
    sorted_total_counts = total_counts.sort_values(ascending=False)[:nb_taxa]  # Sort the top x total counts in descending order
    other_taxa = total_counts.sort_values(ascending=False)[nb_taxa:] # Sort the rest
    
    # Calculate the values for the pie chart
    values = sorted_total_counts.values
    other_percent = other_taxa.sum()
    other_series = pd.Series(other_percent, index=['Others'])
    filtered_families = pd.Series(values, index=sorted_total_counts.index)
    filtered_families = pd.concat([filtered_families, other_series]) # Merge dataframe of top x and 'Others' categories

    # Pie chart of taxonomic composition
    plt.figure(figsize=(12, 8))
    plt.pie(filtered_families, autopct='%1.1f%%', startangle=140) #labels=filtered_families.index to have label above each categorie
    plt.legend(filtered_families.index, title=taxo_level, loc='best', bbox_to_anchor=(1, 1))
    plt.title(f'Proportion of taxonomic category {taxo_level} in all samples')
    plt.axis('equal')
    format_file = input("Which file format would you like to save the plot ? SVG or PDF or PNG")
    if format_file == 'PDF':
        plt.savefig("piechart.pdf", format='pdf', pad_inches=0.2)
    elif format_file == 'SVG':
        plt.savefig("piechart.svg", format='svg', pad_inches=0.2)
    elif format_file == 'PNG':
        plt.savefig("piechart.png", format='png', pad_inches=0.2)
    plt.show()
        
    
def venn_diagram(database): # Function to create a Venn diagram of taxonomic overlap between samples
    # Extract column names of samples
    sequence_index = list(database.columns).index("sequence")
    samples_columns = list(database.columns[sequence_index + 1:]) # Selection of the columns corresponding to the samples (after the “sequence” column)
    print("Available samples : ")
    for sample in samples_columns:
        print(sample)

    sample_1 = input("First sample you want to analyze :")
    sample_2 = input("Second sample you want to analyze :")
    sample_3 = input("Third sample you want to analyze (not mandatory) :")
    
    # New list to store wanted samples
    selected_samples = [sample_1, sample_2]
    if sample_3:
        selected_samples.append(sample_3)
    
    taxo_level = input(f'What taxonomic level do you want ? \n {taxonomy_levels} or ASV \nYour choice : ')
    if taxo_level == 'ASV':
        taxo_level="ASVNumber"
        
    # Loop through each sample column to identify taxa present in each sample
    taxo_samples = [set(database[database[colonne].notna() & (database[colonne] > 0)][taxo_level]) for colonne in selected_samples]
    
    # Choose the appropriate Venn diagram based on the number of samples
    if not sample_3 :
        venn2(taxo_samples, (selected_samples[0], selected_samples[1]))
    else :
        venn3(taxo_samples, (selected_samples[0], selected_samples[1], selected_samples[2]))
    
    if taxo_level == 'ASVNumber':
        plt.title(f'Venn diagram of ASV in samples')
    else:
        plt.title(f'Venn diagram of taxonomic category {taxo_level} in samples')
    plt.tight_layout()
    format_file = input("Which file format would you like to save the plot ? SVG or PDF or PNG")
    if format_file == 'PDF':
        plt.savefig("venn_diagram.pdf", format='pdf', pad_inches=0.5)
    elif format_file == 'SVG':
        plt.savefig("venn_diagram.svg", format='svg', pad_inches=0.5)
    elif format_file == 'PNG':
        plt.savefig("venn_diagram.png", format='png', pad_inches=0.5)
    plt.show()
        
        
def composition_diagram(database):  # Function to create a composition diagram showing the distribution of taxonomic categories in samples
    # Extract column names of samples
    sequence_index = list(database.columns).index("sequence")
    samples_columns = database.columns[sequence_index + 1:]

    taxo_level = input(f'What taxonomic level do you want ? \n {taxonomy_levels} \nYour choice : ')
            
    taxo_samples = {} # Dictionary to store taxonomic categories and their frequencies for each sample
    
    # Iterate through each sample column to calculate the frequency of each taxonomic category
    for col in samples_columns:
        taxo_samples[col] = database[database[col].notna() & (database[col] > 0)][taxo_level].value_counts(normalize=True)
    data = pd.DataFrame(taxo_samples)

    # Bar plot to visualize the distribution of taxonomic categories in samples
    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(111)
    data.plot(kind='bar', ax=ax)
    plt.xlabel(taxo_level)
    plt.ylabel('Frequency in samples')
    plt.title(f'Diagram of {taxo_level} taxonomic category distribution in samples')
    plt.legend(title="Samples", loc='best', bbox_to_anchor=(1, 1))
    plt.xticks(rotation=90, ha='right')
    plt.tight_layout()
    format_file = input("Which file format would you like to save the plot ? SVG or PDF or PNG")
    if format_file == 'PDF':
        plt.savefig("diagram_taxo_composition.pdf", format='pdf', bbox_inches='tight', pad_inches=0.2)
    elif format_file == 'SVG':
        plt.savefig("diagram_taxo_composition.svg", format='svg', bbox_inches='tight', pad_inches=0.2)
    elif format_file == 'PNG':
        plt.savefig("diagram_taxo_composition.png", format='png', bbox_inches='tight', pad_inches=0.2)
    plt.show()
    
    
def network(asv_taxo): # Function to create a network graph based on taxonomic similarity
    taxo_level = input(f'What taxonomic level do you want ? \n {taxonomy_levels} \nYour choice : ')

    # Filter out rows where the taxonomic level is empty (NA)
    filtered_asv_taxo = asv_taxo[(asv_taxo[taxo_level] != 'NA') & (asv_taxo[taxo_level].notna()) & (asv_taxo[taxo_level] != '')]
    G = nx.Graph() # Undirected graph
       
    # Add nodes corresponding to ASV numbers
    for asv in filtered_asv_taxo["ASVNumber"]:
        G.add_node(asv)

    # Add edges between ASVs based on similarity for the specified taxonomic level
    for i, row1 in filtered_asv_taxo.iterrows():
        for j, row2 in filtered_asv_taxo.iterrows():
            if i < j:
                if row1[taxo_level] == row2[taxo_level] and row1["ASVNumber"] != row2["ASVNumber"]:
                    G.add_edge(row1["ASVNumber"], row2["ASVNumber"])
                    
    nx.draw(G, with_labels=True)
    format_file = input("Which file format would you like to save the plot ? SVG or PDF or PNG")
    if format_file == 'PDF':
        plt.savefig("network.pdf", format='pdf', pad_inches=0.2)
    elif format_file == 'SVG':
        plt.savefig("network.svg", format='svg', pad_inches=0.2)
    elif format_file == 'PNG':
        plt.savefig("network.png", format='png', pad_inches=0.2)
    plt.show()
            


####################
####### MAIN #######
####################

if __name__ == "__main__":
    asv_info=sys.argv[1] # Path to "asv.csv" file, output from dada2 [mandatory]
    asv_taxo=sys.argv[2] # Path to "taxo.csv" file, output from dada2 [mandatory]
    database=sys.argv[3] # Path to "database.csv" file, output from dada2 [mandatory]
    db_used=sys.argv[4]
    
    if len(sys.argv) >= 6:
        condition = sys.argv[5] # Path to the file (csv or excel) containing the names of the samples and the conditions to which the samples belong
        extensions_excel = ['.xls', '.xlsx', '.xlsm', '.xlsb', '.odf', '.ods', '.odt']
        if condition.endswith('.csv'):
            condition = pd.read_csv(condition, sep=';')
        elif any(condition.endswith(ext) for ext in extensions_excel):
            condition = pd.read_excel(condition)
        else:
            raise ValueError("File must be in CSV or any Excel file .")
    else:
        condition = None
    
    asv_info = pd.read_csv(asv_info)
    asv_taxo = pd.read_csv(asv_taxo)
    database = pd.read_csv(database)
    
    
    ### Data normalization ###
    scaler_choice = input("Which type would you like for data normalization ? Standardization / Rarefaction / CSS / DESeq2 :  ")
    if scaler_choice == 'Standardization':
        standard_choice = input("Which scaler would you like for standardization ? MinMax (range 0-1) / StandardScaler (mean0-variance1) :  ")
        if standard_choice == 'MinMax':
            scaler = MinMaxScaler()
            
            # Normalization of asv_info
            asv_info_normalized = asv_info.copy()
            columns_numeric = asv_info.columns[1:]
            asv_info_normalized[columns_numeric] = scaler.fit_transform(asv_info_normalized[columns_numeric])
            
            # Normalization of database
            database_normalized = database.copy()
            sequence_index = list(database.columns).index("sequence")
            columns_numeric = database.columns[sequence_index + 1:]
            database_normalized[columns_numeric] = scaler.fit_transform(database_normalized[columns_numeric])

        elif standard_choice == 'StandardScaler':
            scaler = StandardScaler()
            
            # Normalization of asv_info
            asv_info_normalized = asv_info.copy()
            columns_numeric = asv_info.columns[1:]
            asv_info_normalized[columns_numeric] = scaler.fit_transform(asv_info_normalized[columns_numeric])
            
            # Find the smallest value in the normalized DataFrame to add this value to all columns to have only positive values
            min_value = asv_info_normalized[columns_numeric].min().min()
            # Calculate the adjustment value to add (to make the smallest number non-negative)
            if min_value < 0:
                adjustment_value = -min_value
            else:
                adjustment_value = 0
            asv_info_normalized[columns_numeric] = asv_info_normalized[columns_numeric] + adjustment_value # Add this adjustment value to all numeric columns
            
            # Normalization of database
            database_normalized = database.copy()
            sequence_index = list(database.columns).index("sequence")
            columns_numeric = database.columns[sequence_index + 1:]
            database_normalized[columns_numeric] = scaler.fit_transform(database_normalized[columns_numeric])
            
            # Find the smallest value in the normalized DataFrame to add this value to all columns to have only positive values
            min_value = database_normalized[columns_numeric].min().min()
            if min_value < 0:
                adjustment_value = -min_value
            else:
                adjustment_value = 0
            database_normalized[columns_numeric] = database_normalized[columns_numeric] + adjustment_value

    
    
    elif scaler_choice == 'Rarefaction':
        # Normalization of asv_info
        asv_info_normalized = asv_info.copy()
        columns_numeric = asv_info.columns[1:]
        asv_info_normalized[columns_numeric] = preprocessing.normalize(asv_info_normalized[columns_numeric], axis=0)

        # Normalization of database
        database_normalized = database.copy()
        sequence_index = list(database.columns).index("sequence")
        columns_numeric = database.columns[sequence_index + 1:]
        database_normalized[columns_numeric] = preprocessing.normalize(database_normalized[columns_numeric], axis=0)
    
    elif scaler_choice == 'CSS':
        # Normalization of asv_info
        asv_info_normalized = asv_info.copy()
        cumsum_data_asv_info = np.cumsum(asv_info_normalized.iloc[:, 1:], axis=0) # Calculate the cumulative sums for each sample
        medians_asv_info = cumsum_data_asv_info.median() #Calculate the medians of the cumulative sums for each sample
        normalized_data_asv_info = asv_info_normalized.iloc[:, 1:].div(medians_asv_info) # Normalize the data by dividing each count by the corresponding sample's median cumulative sum
        asv_info_normalized = pd.concat([asv_info_normalized.iloc[:, 0], normalized_data_asv_info], axis=1) # Add the ASV names as the first column
        
        # Normalization of database
        database_normalized = database.copy()
        sequence_index = list(database.columns).index("sequence")
        cumsum_data_database = np.cumsum(database_normalized.iloc[:, sequence_index + 1:], axis=0) # Calculate the cumulative sums for each sample
        medians_database = cumsum_data_database.median() #Calculate the medians of the cumulative sums for each sample
        normalized_data_database = database_normalized.iloc[:, sequence_index + 1:].div(medians_database) # Normalize the data by dividing each count by the corresponding sample's median cumulative sum
        database_normalized = pd.concat([database_normalized.iloc[:, :sequence_index + 1], normalized_data_database], axis=1) # Add the other columns
    
    elif scaler_choice == 'DESeq2':
        
        # Normalization of asv_info
        asv_info_normalized = asv_info.copy()
        normalized_data_asv_info = asv_info.iloc[:, 1:]
        normalized_data_asv_info.index = asv_info.iloc[:, 0]
        normalized_data_asv_info_t = normalized_data_asv_info.transpose()
        deseq2_counts_asv_info, size_factors_asv_info = deseq2_norm(normalized_data_asv_info_t) # Use of PyDESeq2 python package
        deseq2_counts_asv_info = deseq2_counts_asv_info.transpose()
        asv_info_normalized = deseq2_counts_asv_info.reset_index()
        
        # Normalization of database
        database_normalized = database.copy()
        sequence_index = list(database.columns).index("sequence")
        normalized_data_database = database.iloc[:, sequence_index + 1:]
        normalized_data_database_t = normalized_data_database.transpose()
        deseq2_counts_database, size_factors_database = deseq2_norm(normalized_data_database_t) # Use of PyDESeq2 python package
        deseq2_counts_database = deseq2_counts_database.transpose()
        database_normalized = pd.concat([database_normalized.iloc[:, :sequence_index + 1], deseq2_counts_database], axis=1) # Add the other columns

    else :
        print("Normalization choice not supported.")
    
    
    if db_used == 'PR2':
        taxonomy_levels = ['Domain', 'Supergroup', 'Division', 'Subdivision', 'Class', 'Order', 'Family', 'Genus', 'Species']
    else :
        taxonomy_levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus']

    
    keep_going = 'yes'
    while keep_going.lower() == 'yes':
        analyse = input("What would you like to study? (alpha_diversity / beta_diversity / graph_taxo_composition):")
        
        if analyse == 'alpha_diversity':
            analyse_alpha = input("What you want to know about alpha diversity ? (diversity_one_sample / diversity_all_samples / statistical_test / alpha_diversity_graph / alpha_diversity_boxplot_condition) : ")
            if analyse_alpha == 'diversity_one_sample':
                alpha_diversity_one(asv_info)
            elif analyse_alpha == 'diversity_all_samples':
                alpha_diversity_all(asv_info)
            elif analyse_alpha == 'statistical_test':
                statistical_test_alpha(asv_info, condition)
            elif analyse_alpha == 'alpha_diversity_graph':
                alpha_graph(asv_info)
            elif analyse_alpha == 'alpha_diversity_boxplot_condition':
                alpha_graph_condition(asv_info, condition)
            else:
                print("Analysis not supported.")
                exit()
            
        elif analyse == 'beta_diversity':
            analyse_beta = input("What you want to know about beta diversity ? (beta_diversity / statistical_test / beta_diversity_graph) : ")
            if analyse_beta == 'beta_diversity':
                beta_diversity = beta_diversity_function(asv_info_normalized)
                beta_diversity_all(beta_diversity)
            elif analyse_beta == 'statistical_test':
                beta_diversity = beta_diversity_function(asv_info_normalized)
                statistical_test_beta(beta_diversity, asv_info_normalized, condition)
            elif analyse_beta == 'beta_diversity_graph':
                beta_diversity = beta_diversity_function(asv_info_normalized)
                beta_diversity_graph(beta_diversity, condition, asv_info_normalized)
            else:
                print("Analysis not supported.")
                exit()
            
        elif analyse == 'graph_taxo_composition':
            composition_level = input("What kind of representation would you like ? (barplot_one_sample / barplot_all_samples / barplot_all_samples_condition / barplot_all_samples_two_conditions / heatmap / piechart / venn_diagram / diagram_taxo_composition /network ): ")
            if composition_level == 'barplot_one_sample':
                barplot_one_sample(database)
            if composition_level == 'barplot_all_samples':
                barplot_all_samples(database)
            if composition_level == 'barplot_all_samples_condition':
                barplot_all_samples_condition(database, condition, asv_info_normalized)
            if composition_level == 'barplot_all_samples_two_conditions':
                barplot_all_samples_two_conditions(database, condition, asv_info_normalized)
            if composition_level == 'heatmap':
                heatmap(database)
            if composition_level == 'piechart':
                piechart(database)
            if composition_level == 'venn_diagram':
                venn_diagram(database_normalized)
            if composition_level == 'diagram_taxo_composition':
                composition_diagram(database)
            if composition_level == 'network':
                network(asv_taxo)

        else:
            print("Analysis not supported.")
            
        keep_going = input("Do you want to run another analysis ? (yes / no) ")

