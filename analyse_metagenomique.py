import sys
import ast
import random
import pandas as pd
import skbio
import seaborn as sns
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import skbio.diversity.alpha as alpha
from skbio import DistanceMatrix
from skbio.tree import nj
from skbio.diversity import beta_diversity
from skbio.stats.ordination import pcoa
from skbio.stats.distance import permanova
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.manifold import MDS
from matplotlib import cm
from matplotlib_venn import venn2, venn3
#from matplotlib.patches import Ellipse
from scipy.stats import f_oneway, tukey_hsd
from statsmodels.stats.multicomp import pairwise_tukeyhsd



def one_condition_struct(asv_info, condition): # Function to create a structure containing different conditions with associated samples
    # Check if a condition file is provided
    if condition is None:
        print("Error: Please provide a condition file.")
        return

    print("Name of conditions available :")
    available_condition = list(condition.columns)
    for name in available_condition:
        print(name)
    column_sample_name = input("What is the exact name of the column with the sample names? : OTUNumber ")
    column_condition_name = input("What is the name of the column with the names of the conditions to which the samples belong? : Condition ")
    
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
    column_sample_name = input("What is the exact name of the column with the sample names? : OTUNumber ")
    column_condition_first_name = input("What is the name of the column with the names of the conditions to which the samples belong? : Condition ")
    column_condition_second_name = input("What is the name of the column with the names of the conditions to which the samples belong? : Condition2 ")
    
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

    return conditions


#####################
## Alpha Diversity ##
#####################


def alpha_diversity_one(asv_info): # Function to calculate alpha diversity for a given sample
    print("Available samples :")
    for column in asv_info.columns[1:]:
        print(column)
    sample_alpha = input("Which sample do you want alpha diversity ? : ")
    counts = asv_info[sample_alpha]     # Retrieve counts for the selected sample
    alpha_index = input("Which alpha diversity index do you want to calculate ? (shannon / simpson / chao / observed): ")
    if alpha_index == 'shannon':
        alpha_diversity = skbio.diversity.alpha.shannon(counts)
        print("-- Shannon Alpha Diversity for sample ", sample_alpha, " : ", alpha_diversity)
        
    elif alpha_index == 'simpson':
        alpha_diversity = skbio.diversity.alpha.simpson(counts)
        print("-- Simpson Alpha Diversity for sample ", sample_alpha, " : ", alpha_diversity)
        
    elif alpha_index == 'chao':
        alpha_diversity = skbio.diversity.alpha.chao1(counts)
        print("-- Chao Alpha Diversity for sample ", sample_alpha, " : ", alpha_diversity)
        
    elif alpha_index == 'observed':
        asv_sample = asv_info[['ASVNumber', sample_alpha]].loc[asv_info[sample_alpha] > 0]
        asv_names = asv_sample['ASVNumber'].tolist()
        asv_number = len(asv_names)
        print("-- Observed Alpha Diversity for sample ",sample_alpha, " : ",  asv_number)
    
    else:
        print("Alpha diversity index not supported.")
        exit()
    
   
def alpha_diversity_all(asv_info): # Function to calculate alpha diversity for all samples
    alpha_diversity_all = {}
    alpha_index = input("Which alpha diversity index do you want to calculate ? (shannon / simpson / chao / observed): ")
    if alpha_index == 'shannon':
        print("-- Shannon Alpha diversity --")
        for column in asv_info.columns[1:]:  # Ignore the first column (ASVNumber) to have samples columns
            counts = asv_info[column] # Retrieve counts for the selected sample
            alpha_diversity_all[column] = {}
            alpha_diversity_all[column]['shannon'] = skbio.diversity.alpha.shannon(counts)
        for column, diversity in alpha_diversity_all.items():
            print(column, " : ", diversity['shannon'])
            
    elif alpha_index == 'simpson':
        print("-- Simpson Alpha diversity --")
        for column in asv_info.columns[1:]:
            counts = asv_info[column]
            alpha_diversity_all[column] = {}
            alpha_diversity_all[column]['simpson'] = skbio.diversity.alpha.simpson(counts)
        for column, diversity in alpha_diversity_all.items():
            print(column, " : ", diversity['simpson'])
            
    elif alpha_index == 'chao':
        print("-- Chao Alpha diversity --")
        for column in asv_info.columns[1:]:
            counts = asv_info[column]
            alpha_diversity_all[column] = {}
            alpha_diversity_all[column]['chao'] = skbio.diversity.alpha.chao1(counts)
        for column, diversity in alpha_diversity_all.items():
            print(column, " : ", diversity['chao'])
            
    elif alpha_index == 'observed':
        print("-- Observed Alpha diversity --")
        asv_all_samples = {}
        for column in asv_info.columns[1:]:
            asv_samples = asv_info[['ASVNumber', column]].loc[asv_info[column] > 0] # Keep only rows of the df where the value in the current column is greater than 0
            asv_names = asv_samples['ASVNumber'].tolist()
            asv_all_samples[column] = len(asv_names) # Number of observed ASVs for the current sample
        for sample, asv in asv_all_samples.items():
            print(sample, " : ", asv)
    else:
        print("Alpha diversity index not supported.")
        exit()
    

def statistical_test_alpha(asv_info, condition): # Function to perform statistical tests on alpha diversity
    conditions=one_condition_struct(asv_info, condition)
    sample_counts = {} # Dictionary to store counts of ASVs for each sample
    sample_data = {} # Dictionary to store ASV data for each condition
    
    for condition,samples_list in conditions.items(): # Looping through each condition and its associated samples
        for sample in samples_list:
            if sample not in sample_counts:
                sample_counts[sample] = []
            sample_counts[sample].extend(asv_info[sample]) # Extending the ASV counts for each sample
        
        for sample_cond in conditions[condition]:
            if condition not in sample_data:
                sample_data[condition] = []
            sample_data[condition].extend(sample_counts[sample_cond]) # Extracting ASV data for each condition based on the associated samples
        
    # Creating a list to hold ASV data for each condition
    condition_data = [[] for j in range(len(sample_data))]
    for i, (cond, data) in enumerate(sample_data.items()):
        condition_data[i] = data
    print("sample_counts")
    print(sample_counts)
    print("sample_data")
    print(sample_data)
    f_statistic, p_value = f_oneway(*condition_data) # One-way ANOVA test
    print("ANOVA test result between all conditions : ")
    print("Statistic F : ", f_statistic)
    print("p-value : ", p_value)
        
    post_hoc = input("Do you want to make post-hoc Tuckey test ? Yes/No : ") # 95% test
    if post_hoc == 'Yes':
        condition_data_flat = np.concatenate(condition_data) # Flattening the ASV data
        conditions_list = [] # Creating a list of conditions for Tukey test
        
        for cond, data in sample_data.items():
            for value in data:
                conditions_list.append(cond) # Add conditions for each ASV

        tukey_result = pairwise_tukeyhsd(endog=condition_data_flat, groups=conditions_list)
        print("Result of Tukey test : ")
        print(tukey_result)

def alpha_graph(asv_info): # Function to create alpha diversity scatter plot
    alpha_index = input("Which alpha diversity index would you like ? (shannon / simpson / chao / observed) : ")
    alpha_diversity = {} # Dictionary to store alpha diversity values for each sample and index
    for column in asv_info.columns[1:]:
        counts = asv_info[column] # Extract counts for the sample
        if alpha_index == 'shannon':
            alpha_diversity.setdefault('shannon', {})[column] = skbio.diversity.alpha.shannon(counts)
        elif alpha_index == 'simpson':
            alpha_diversity.setdefault('simpson', {})[column] = skbio.diversity.alpha.simpson(counts)
        elif alpha_index == 'chao':
            alpha_diversity.setdefault('chao', {})[column] = skbio.diversity.alpha.chao1(counts)
        elif alpha_index == 'observed':
            asv_samples = asv_info[['ASVNumber', column]].loc[asv_info[column] > 0] # Keep only rows of the df where the value in the current column is greater than 0
            asv_names = asv_samples['ASVNumber'].tolist()
            alpha_diversity.setdefault('observed', {})[column] = len(asv_names) # Number of observed ASVs for the current sample
        else:
            print("Alpha diversity index not supported.")
            exit()

    # Create the figure
    plt.figure(figsize=(12, 10))
    plt.scatter(alpha_diversity[alpha_index].keys(), alpha_diversity[alpha_index].values(), color='blue')
    plt.xlabel('Samples')
    plt.ylabel(f'{alpha_index.capitalize()} Alpha Diversity')
    plt.title(f'{alpha_index.capitalize()} Alpha Diversity according to samples')
    plt.xticks(rotation=90)
    plt.tight_layout()
    format_file = input("Which file format would you like to save the plot ? SVG or PDF ")
    if format_file == 'PDF':
        plt.savefig("alpha_diversity.pdf", format='pdf', pad_inches=0.2)
    elif format_file == 'SVG':
        plt.savefig("alpha_diversity.svg", format='svg', pad_inches=0.2)
    plt.show()


def alpha_graph_condition(asv_info, condition): # Function to create alpha diversity boxplots grouped by conditions
    alpha_index = input("Which alpha diversity index would you like? (shannon / simpson / chao / observed): ")
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
            elif alpha_index == 'chao':
                alpha_results[cond].append(skbio.diversity.alpha.chao1(counts))
            elif alpha_index == 'observed':
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
    bp = plt.boxplot(alpha_results.values(), labels=alpha_results.keys(), patch_artist=True, capprops={'linewidth': 0.0}) #, whiskerprops={'linewidth': 0.0} pour enlever moustache
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
    format_file = input("Which file format would you like to save the plot ? SVG or PDF ")
    if format_file == 'PDF':
        plt.savefig("boxplot_alpha_div.pdf", format='pdf', pad_inches=0.2)
    elif format_file == 'SVG':
        plt.savefig("boxplot_alpha_div.svg", format='svg', pad_inches=0.2)
    plt.show()



####################
## Beta Diversity ##
####################


def beta_diversity_function(asv_info): # Function to calculate beta diversity
    count_table = {}  # Dictionary to store count tables for each sample
    for column in asv_info.columns[1:]:
        asv_counts = {} # Dictionary to store counts for each ASV in the sample
        for index, row in asv_info.iterrows():
            count = row[column] # Get the count of this ASV in this sample
            asv_counts[index] = count
        count_table[column] = asv_counts
    df_counts = pd.DataFrame(count_table).T.fillna(0)  # Transpose so that samples are rows
        
    beta_index = input("Which alpha diversity indices would you like ? : (euclidean/cityblock/braycurtis/canberra/chebyshev/correlation/cosine/dice/hamming/jaccard/kulsinski/matching/minkowski/rogerstanimoto/russellrao/seuclidean/sokalmichener/sokalsneath/sqeuclidean/yule)")
    if beta_index not in ['braycurtis', 'canberra', 'chebyshev', 'cityblock', 'correlation', 'cosine', 'dice', 'euclidean', 'hamming', 'jaccard', 'kulsinski', 'matching', 'minkowski', 'rogerstanimoto', 'russellrao', 'seuclidean', 'sokalmichener', 'sokalsneath', 'sqeuclidean', 'yule']:
        print("Beta diversity index not supported. The available index are : braycurtis, canberra, chebyshev, cityblock, correlation, cosine, dice, euclidean, hamming, jaccard, kulsinski, matching, minkowski, rogerstanimoto, russellrao, seuclidean, sokalmichener, sokalsneath, sqeuclidean, yule ")
        exit()
    
    # Calculate beta diversity using the selected index and count table
    beta_diversity = skbio.diversity.beta_diversity(beta_index, df_counts)
    sample_names = list(asv_info.columns[1:]) # Extract sample name
    beta_diversity.index = sample_names # Set the IDs of beta diversity results samples to the sample names
    
    return df_counts, beta_diversity
    
    
def beta_diversity_all(beta_diversity): # Function to display beta diversity for all samples
    print("-- Beta diversity :")
    print(beta_diversity)

        
        
def statistical_test_beta(df_counts, beta_diversity): # Function to perform a statistical test (permanova) on beta diversity
    # Extract groups from sample names
    groups = [nom.split('_')[0] for nom in df_counts.index]
    permanova_result = permanova(distance_matrix=beta_diversity, grouping=groups)
    print("-- Permanova test results")
    print("Permanova statistic : ", permanova_result['test statistic'])
    print("p-value : ", permanova_result['p-value'])
    print("Number of permutations : ", permanova_result['number of permutations'])
    
    
def beta_diversity_graph(beta_diversity, condition, asv_info): # Function to create beta diversity visualizations
    beta_representation = input("Which beta diversity representation would you like ? (heatmap / NMDS / NMDS_condition / PCoA / PCoA_condition): ")
    if beta_representation == 'heatmap': # Create a heatmap representation of beta diversity
        plt.figure(figsize=(10, 8))
        sns.heatmap(beta_diversity.to_data_frame(), cmap="rainbow_r", annot=True, fmt=".2f", linewidths=.5)
        plt.title("Beta diversity heatmap")
        plt.xlabel("Samples")
        plt.ylabel("Samples")
        format_file = input("Which file format would you like to save the plot ? SVG or PDF ")
        if format_file == 'PDF':
            plt.savefig("heatmap_beta_diversity.pdf", format='pdf', pad_inches=0.2)
        elif format_file == 'SVG':
            plt.savefig("heatmap_beta_diversity.svg", format='svg', pad_inches=0.2)
        plt.show()
        
    elif beta_representation == 'NMDS':
        mds = MDS(metric=False, random_state=0) # MDS object which perform nonmetric MDS with reproducibility
        beta_diversity_array = np.array(beta_diversity[:][:])
        mds_results = mds.fit_transform(beta_diversity_array) # Performs MDS transformation
        stress = mds.stress_
        print(stress)
        if stress <= 0.2:
            colors = plt.cm.rainbow(np.linspace(0, 1, len(beta_diversity_array))) # Generates colors for each sample
            
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
            plt.legend(handles=legend_handles, title='Echantillons', loc='best', bbox_to_anchor=(1, 1))
            plt.tight_layout()
            plt.xlabel('Axis 1')
            plt.ylabel('Axis 2')
            plt.title('NMDS Plot')
            format_file = input("Which file format would you like to save the plot ? SVG or PDF ")
            if format_file == 'PDF':
                plt.savefig("NMDS.pdf", format='pdf', bbox_inches='tight', pad_inches=0.2)
            elif format_file == 'SVG':
                plt.savefig("NMDS.svg", format='svg', bbox_inches='tight', pad_inches=0.2)
            plt.show()
        else :
            print("The stress variable is ", stress, ". It's greater than 0.2. Perform a PCoA analysis instead.")
            pcoa_results = pcoa(beta_diversity) # Perform PCoA
            sample_names = list(asv_info.columns[1:]) # Extract sample name
            pcoa_results.samples.index = sample_names # Set the IDs of PCoA results samples to the sample names
            colors = plt.cm.rainbow(np.linspace(0, 1, len(sample_names))) # Generate a list of colors
            legend_handles = []
            plt.figure(figsize=(12, 8))
            
            for sample, color in zip(sample_names, colors): # Iterate over sample names and corresponding colors
                x = pcoa_results.samples.loc[sample, 'PC1']
                y = pcoa_results.samples.loc[sample, 'PC2']
                handle = plt.scatter(x, y, color=color, label=sample)
                legend_handles.append(handle)
                
            plt.legend(handles=legend_handles, title='Echantillons', loc='best', bbox_to_anchor=(1, 1))
            plt.tight_layout()
            plt.title('PCoA Plot')
            plt.xlabel('Axis 1')
            plt.ylabel('Axis 2')
            format_file = input("Which file format would you like to save the plot ? SVG or PDF ")
            if format_file == 'PDF':
                plt.savefig("PCoA.pdf", format='pdf', bbox_inches='tight', pad_inches=0.2)
            elif format_file == 'SVG':
                plt.savefig("PCoA.svg", format='svg', bbox_inches='tight', pad_inches=0.2)
            plt.show()
        
    elif beta_representation == 'NMDS_condition':
        mds = MDS(metric=False, random_state=0) # MDS object which perform nonmetric MDS with reproducibility
        beta_diversity_array = np.array(beta_diversity[:][:])
        mds_results = mds.fit_transform(beta_diversity_array) # Performs MDS transformation
        stress = mds.stress_
        print(stress)
        if stress <= 0.2:
            mds_results_df = pd.DataFrame(mds_results, columns=['Dimension 1', 'Dimension 2'])
            # Replace row and column numbers with sample names
            sample_names = list(asv_info.columns[1:])
            mds_results_df.index = sample_names
            print(mds_results_df)
            
            # Retrieve the structure containing sample data based on conditions
            conditions=one_condition_struct(asv_info, condition)
            # Extract legend labels from conditions
            legend_labels = []
            for condition,samples_list in conditions.items():
                legend_labels.append(condition)
            
            # Generate random unique colors for each condition
            colors = {}
            for condition in legend_labels:
                color = "#{:06x}".format(random.randint(0, 0xFFFFFF))
                colors[condition] = color
                
            # Assign a color to each sample based on its condition
            sample_colors = [colors[condition] for sample in mds_results_df.index for condition, samples_list in conditions.items() if sample in samples_list]
                  
            # Scatter plot with samples colored by condition
            plt.figure(figsize=(12, 8))
            for condition, color in colors.items():
                indices = [i for i, sample in enumerate(mds_results_df.index) if sample in conditions[condition]]
                plt.scatter(mds_results_df.iloc[indices, 0], mds_results_df.iloc[indices, 1], c=color, label=condition)
            plt.legend(title='Conditions', loc='best', bbox_to_anchor=(1, 1))
            plt.xlabel('Axis 1')
            plt.ylabel('Axis 2')
            plt.title('NMDS Plot')
            format_file = input("Which file format would you like to save the plot ? SVG or PDF ")
            if format_file == 'PDF':
                plt.savefig("NMDS_condition.pdf", format='pdf', bbox_inches='tight', pad_inches=0.2)
            elif format_file == 'SVG':
                plt.savefig("NMDS_condition.svg", format='svg', bbox_inches='tight', pad_inches=0.2)
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
            
            # Generate random unique colors for each condition
            colors = {}
            for condition in legend_labels:
                color = "#{:06x}".format(random.randint(0, 0xFFFFFF))
                colors[condition] = color
                
            # Assign a color to each sample based on its condition
            sample_colors = [colors[condition] for sample in pcoa_results.samples.index for condition, samples_list in conditions.items() if sample in samples_list]
                  
            # Scatter plot with samples colored by condition
            plt.figure(figsize=(12, 8))
            for condition, color in colors.items():
                indices = [i for i, sample in enumerate(pcoa_results.samples.index) if sample in conditions[condition]] # Find indices of samples belonging to the current condition
                plt.scatter(pcoa_results.samples.iloc[indices, 0], pcoa_results.samples.iloc[indices, 1], c=color, label=condition)  # Scatter plot the samples belonging to the current condition using PC1 and PC2 coordinates
                
            plt.legend(title='Conditions', loc='best', bbox_to_anchor=(1, 1))
            plt.tight_layout()
            plt.title('PCoA Plot')
            plt.xlabel('Axis 1')
            plt.ylabel('Axis 2')
            format_file = input("Which file format would you like to save the plot ? SVG or PDF ")
            if format_file == 'PDF':
                plt.savefig("PCoA_condition.pdf", format='pdf', bbox_inches='tight', pad_inches=0.2)
            elif format_file == 'SVG':
                plt.savefig("PCoA_condition.svg", format='svg', bbox_inches='tight', pad_inches=0.2)
            plt.show()
        
        
    elif beta_representation == 'PCoA': # Create a PCoA plot of beta diversity, colored by samples
        pcoa_results = pcoa(beta_diversity) # Perform PCoA
        sample_names = list(asv_info.columns[1:]) # Extract sample name
        pcoa_results.samples.index = sample_names # Set the IDs of PCoA results samples to the sample names
        colors = plt.cm.rainbow(np.linspace(0, 1, len(sample_names))) # Generate a list of colors
        legend_handles = []
        plt.figure(figsize=(12, 8))
        
        for sample, color in zip(sample_names, colors): # Iterate over sample names and corresponding colors
            x = pcoa_results.samples.loc[sample, 'PC1']
            y = pcoa_results.samples.loc[sample, 'PC2']
            handle = plt.scatter(x, y, color=color, label=sample)
            legend_handles.append(handle)
            
        plt.legend(handles=legend_handles, title='Echantillons', loc='best', bbox_to_anchor=(1, 1))
        plt.tight_layout()
        plt.title('PCoA Plot')
        plt.xlabel('Axis 1')
        plt.ylabel('Axis 2')
        format_file = input("Which file format would you like to save the plot ? SVG or PDF ")
        if format_file == 'PDF':
            plt.savefig("PCoA.pdf", format='pdf', bbox_inches='tight', pad_inches=0.2)
        elif format_file == 'SVG':
            plt.savefig("PCoA.svg", format='svg', bbox_inches='tight', pad_inches=0.2)
        plt.show()
    
    elif beta_representation == 'PCoA_condition': # Create a PCoA plot of beta diversity grouped by conditions
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
        
        # Generate random unique colors for each condition
        colors = {}
        for condition in legend_labels:
            color = "#{:06x}".format(random.randint(0, 0xFFFFFF))
            colors[condition] = color
            
        # Assign a color to each sample based on its condition
        sample_colors = [colors[condition] for sample in pcoa_results.samples.index for condition, samples_list in conditions.items() if sample in samples_list]
              
        # Scatter plot with samples colored by condition
        plt.figure(figsize=(12, 8))
        for condition, color in colors.items():
            indices = [i for i, sample in enumerate(pcoa_results.samples.index) if sample in conditions[condition]] # Find indices of samples belonging to the current condition
            plt.scatter(pcoa_results.samples.iloc[indices, 0], pcoa_results.samples.iloc[indices, 1], c=color, label=condition)  # Scatter plot the samples belonging to the current condition using PC1 and PC2 coordinates
            
            #-------------------------------#
            # Calculate mean and covariance matrix for the points of each condition
            #mean = pcoa_results.samples.iloc[indices, :].mean(axis=0)
            #cov_matrix = pcoa_results.samples.iloc[indices, :].cov()
            
            # Calculate the eigenvalues and eigenvectors of the covariance matrix
            #eigenvalues, eigenvectors = np.linalg.eigh(cov_matrix)
            #print(eigenvalues)
            
            
            #eigenvalues = np.sqrt(eigenvalues)
            #valid_eigenvalues = [eigval for eigval in eigenvalues if eigval > 1e-6] #permet d'enlever les valeurs nulles ou négatives pour faire la racine
            #eigenvalues = np.sqrt(valid_eigenvalues)
            #print(eigenvalues)


            # Choose the angle of the ellipse (in radians) based on the principal component with the largest eigenvalue
            #angle = np.degrees(np.arctan2(eigenvectors[1, 0], eigenvectors[0, 0]))  # Use the first eigenvector
            #angle=np.rad2deg(np.arccos(eigenvectors[0, 0]))

            # Plot the ellipse
            #ellipse = Ellipse(xy=mean, width=eigenvalues[0]*2, height=eigenvalues[1]*2, angle=angle, color=color, alpha=0.2)
            #plt.gca().add_patch(ellipse)
            #---------------------------------#
        
        plt.legend(title='Conditions', loc='best', bbox_to_anchor=(1, 1))
        plt.tight_layout()
        plt.title('PCoA Plot')
        plt.xlabel('Axis 1')
        plt.ylabel('Axis 2')
        format_file = input("Which file format would you like to save the plot ? SVG or PDF ")
        if format_file == 'PDF':
            plt.savefig("PCoA_condition.pdf", format='pdf', bbox_inches='tight', pad_inches=0.2)
        elif format_file == 'SVG':
            plt.savefig("PCoA_condition.svg", format='svg', bbox_inches='tight', pad_inches=0.2)
        plt.show()
        
    else:
        print("Beta diversity representation not supported.")
        exit()
    


#################################
## Graph taxonomic composition ##
#################################


def barplot_one_sample(database): # Function to create a bar plot of taxonomic composition for a given sample
    print("Available samples :")
    available_samples = list(database.columns[database.columns.get_loc('sequence')+1:])
    for sample in available_samples:
        print(sample)
    sample = input("Which sample would you like to analyze ? ")
    taxo_level = input("What taxonomic level do you want ? \nFor assignment with PR2 : Domain;Supergroup;Division;Subdivision;Class;Order;Family;Genus;Species \nFor assignment with Silva : Kingdom;Phylum;Class;Order;Family;Genus \nYour choice : ")

    taxo_count = database[[taxo_level] + [sample]] # Extract taxonomic counts for the selected sample and taxonomic level
    taxo_sum = taxo_count.groupby(taxo_level).sum() # Group by taxonomic level and sum counts
    taxo_sum[sample] = taxo_sum[sample].astype(float) # Convert count values to float
    taxo_proportions = taxo_sum / taxo_sum.sum() # Calculate proportions of each taxonomic category

    top_20_taxa = taxo_proportions.apply(lambda x: x.nlargest(20), axis=0)
    other_taxa = 1 - top_20_taxa.sum(axis=0) # Calculation of the sum of the proportions of the remaining categories for each sample
    top_20_taxa.loc['Others'] = other_taxa.values # Addition of the 'Others' line with the proportions of the 'Others' category for each sample
    
    # Merge taxonomic categories with proportions less than 0.001 into an "Others" category
    #taxo_others = taxo_proportions[taxo_proportions[sample] < 0.001]
    #proportion_others = taxo_others.sum()
    #taxo_proportions_without_others = taxo_proportions[taxo_proportions[sample] >= 0.001]
    #taxo_proportions_finals = pd.concat([taxo_proportions_without_others, pd.DataFrame([proportion_others], index=['Autres'], columns=[sample])])

    # Bar plot of taxonomic proportions
    plt.figure(figsize=(10, 6))
    #taxo_proportions_finals[sample].sort_values(ascending=False).plot(kind='bar')
    top_20_taxa[sample].sort_values(ascending=False).plot(kind='bar')
    plt.xlabel(taxo_level)
    plt.ylabel('Proportion')
    plt.title(f'Proportion of {taxo_level} taxonomic category in sample {sample}')
    plt.xticks(rotation=90, ha='right')
    plt.tight_layout()
    format_file = input("Which file format would you like to save the plot ? SVG or PDF ")
    if format_file == 'PDF':
        plt.savefig("barplot_one_sample.pdf", format='pdf', pad_inches=0.2)
    elif format_file == 'SVG':
        plt.savefig("barplot_one_sample.svg", format='svg', pad_inches=0.2)
    plt.show()

    
def barplot_all_samples(database): # Function to create a bar plot of taxonomic composition for all samples
    # Select columns corresponding to samples (after the "sequence" column)
    sequence_index = list(database.columns).index("sequence")
    samples_columns = database.columns[sequence_index + 1:]
    taxo_level = input("What taxonomic level do you want ? \nFor assignment with PR2 : Domain;Supergroup;Division;Subdivision;Class;Order;Family;Genus;Species \nFor assignment with Silva : Kingdom;Phylum;Class;Order;Family;Genus \nYour choice : ")
    
    taxo_count = database[[taxo_level] + list(samples_columns)] # Extract taxonomic counts for the selected taxonomic level and all samples
    taxo_sum = taxo_count.groupby(taxo_level).sum() # Group by taxonomic level and sum counts
    taxo_proportions = taxo_sum.div(taxo_sum.sum(axis=0), axis=1) # Calculate proportions of each taxonomic category for each sample
    
    # Selection of the 20 most represented categories for each sample :
    top_20_taxa = taxo_proportions.apply(lambda x: x.nlargest(20), axis=0)
    other_taxa = 1 - top_20_taxa.sum(axis=0) # Calculation of the sum of the proportions of the remaining categories for each sample
    top_20_taxa.loc['Others'] = other_taxa.values # Addition of the 'Others' line with the proportions of the 'Others' category for each sample
    top_20_taxa_transpose = top_20_taxa.transpose() # Transposition to have samples in rows and categories in columns
        
    # Bar plot of taxonomic proportions for all samples
    fig, ax = plt.subplots(figsize=(18, 12))
    top_20_taxa_transpose.plot(kind='bar', stacked=True, colormap='rainbow', ax=ax)
    #taxo_proportions_transpose.plot(kind='bar', stacked=True, colormap='rainbow', ax=ax)
    plt.xlabel('Samples')
    plt.ylabel('Abundance')
    plt.title(f'Abundance of taxonomic category {taxo_level} for all samples')
    plt.legend(title=taxo_level, loc='best', bbox_to_anchor=(1, 1))
    plt.xticks(rotation=90, ha='right')
    plt.tight_layout()
    format_file = input("Which file format would you like to save the plot ? SVG or PDF ")
    if format_file == 'PDF':
        plt.savefig("barplot_composition_all_samples.pdf", format='pdf', bbox_inches='tight', pad_inches=0.2)
    elif format_file == 'SVG':
        plt.savefig("barplot_composition_all_samples.svg", format='svg', bbox_inches='tight', pad_inches=0.2)
    plt.show()
    exit()
        
        

def barplot_all_samples_condition(database, condition, asv_info):
    # Extract column names of samples
    sequence_index = list(database.columns).index("sequence")
    samples_columns = database.columns[sequence_index + 1:]
    
    taxo_level = input("What taxonomic level do you want ? \nFor assignment with PR2 : Domain;Supergroup;Division;Subdivision;Class;Order;Family;Genus;Species \nFor assignment with Silva : Kingdom;Phylum;Class;Order;Family;Genus \nYour choice : ")
    
    conditions=one_condition_struct(asv_info, condition) # Get the list of samples associated with the condition
    
    plot_level = input("Do you want a merged chart with all conditions or separate charts for each condition ? plots / merged_plots / separate_plots ")
    
    if plot_level == 'plots':
        # Initialize a dictionary to store summed abundances for each condition
        condition_abundances = {cond: 0 for cond in conditions}
        
        # Iterate over each condition
        for i, (condition, samples_list) in enumerate(conditions.items()):
            # Select samples corresponding to the condition
            taxo_count = database[[taxo_level] + list(samples_list)]
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

        color_map = sns.color_palette("tab20", len(all_taxons))
        
        for i, (condition_name, taxon_abundances) in enumerate(condition_abundances.items()):
            top_20_taxa = taxon_abundances.sort_values(ascending=False)[:20]
            other_taxa = 1 - top_20_taxa.sum(axis=0)
            top_20_taxa.loc['Others'] = other_taxa # Addition of the 'Others' line with the proportions of the 'Others' category for each condition
            
            # Initialize y_offset to keep track of the vertical position of the next bar
            y_offset = 0
            for j, (taxon, abundance) in enumerate(top_20_taxa.items()):
                if taxon in top_20_taxa:
                    if taxon != 'Others':
                        color = color_map[j] # Assign color based on taxon
                    else:
                        color = 'yellow' # Assign color for 'Others' categorie
                    plt.bar(condition_name, abundance, color=color, bottom=y_offset)  # Plot each taxon with assigned color
                    y_offset += abundance
                    if i == 0:  # Only add handles from the first condition to avoid duplicates
                        legend_handles.append(taxon)
            
        plt.xlabel('Condition')
        plt.ylabel('Total Abundance')
        plt.title('Total Abundance by Condition')
        plt.xticks(rotation=90)
        plt.legend(legend_handles, title=taxo_level, loc='best', bbox_to_anchor=(1, 1))
        plt.tight_layout()
        
        format_file = input("Which file format would you like to save the plot ? SVG or PDF ")
        if format_file == 'PDF':
            plt.savefig("barplot_total_abundance_by_condition.pdf", format='pdf', bbox_inches='tight', pad_inches=0.2)
        elif format_file == 'SVG':
            plt.savefig("barplot_total_abundance_by_condition.svg", format='svg', bbox_inches='tight', pad_inches=0.2)
        plt.show()

    if plot_level == 'merged':
        fig, axs = plt.subplots(1, len(conditions), figsize=(18, 12), sharey=True)
        num_conditions = len(conditions)
        num_cols = len(conditions) #2  # Nombre de colonnes dans la grille
        num_rows = (num_conditions + num_cols - 1) // num_cols  # Calcul du nombre de lignes nécessaire
        
        taxon_colors = {} # Dictionary to store unique colors for each taxon
        
        # Barplots for each condition
        for i, (condition, samples_list) in enumerate(conditions.items()):
            row = i // num_cols  # Ligne actuelle dans la grille
            col = i % num_cols   # Colonne actuelle dans la grille
            
            taxo_count = database[[taxo_level] + list(samples_list)]
            taxo_sum = taxo_count.groupby(taxo_level).sum()
            taxo_proportions = taxo_sum.div(taxo_sum.sum(axis=0), axis=1)
            
            # Select top 20 taxonomic categories
            top_20_taxa = taxo_proportions.sum(axis=1).nlargest(20)
            taxo_proportions = taxo_proportions.transpose() # Transposition pour avoir les échantillons en lignes et les catégories en colonnes
            top_20_taxa = taxo_proportions[top_20_taxa.index]
            other_taxa = 1 - top_20_taxa.sum(axis=1) # Calcul de la somme des proportions des catégories restantes pour chaque échantillon
            top_20_taxa['Others'] = other_taxa.values # Ajout de la ligne 'Others' avec les proportions de la catégorie 'Others' pour chaque échantillon
                        
            # Associate one unique color to each taxa
            for taxon in top_20_taxa.columns:
                if taxon not in taxon_colors:
                # Générer une couleur aléatoire unique pour chaque taxon
                    color = "#{:06x}".format(random.randint(0, 0xFFFFFF))
                    while color in taxon_colors.values():  # Vérifie si la couleur est déjà utilisée
                        color = "#{:06x}".format(random.randint(0, 0xFFFFFF))  # Génère une nouvelle couleur si nécessaire
                    taxon_colors[taxon] = color  # Associer le taxon à la couleur générée
                    
                    #taxon_colors[taxon] = plt.cm.rainbow(len(taxon_colors) / 25)
            
            # Plot the barplot for the current condition
            ax = axs[i]
            top_20_taxa.plot(kind='bar', stacked=True, color=[taxon_colors[taxon] for taxon in top_20_taxa.columns], ax=ax, label=condition, legend=None)
            plt.xlabel('Samples')
            plt.ylabel('Abundance')
            plt.title(condition)
            plt.xticks(rotation=90, ha='right')
        
        plt.tight_layout()
        # Créer une légende commune
        legend_elements = [plt.Line2D([0], [0], color=taxon_colors[taxon], lw=4, label=taxon) for taxon in taxon_colors]
        plt.legend(handles=legend_elements, title=taxo_level, loc='best', bbox_to_anchor=(1, 1))

        # Sauvegarder et afficher le plot
        format_file = input("Which file format would you like to save the plot ? SVG or PDF ")
        if format_file == 'PDF':
            plt.savefig("barplots_composition.pdf", format='pdf', bbox_inches='tight', pad_inches=0.2)
        elif format_file == 'SVG':
            plt.savefig("barplots_composition.svg", format='svg', bbox_inches='tight', pad_inches=0.2)
        plt.show()


    
    if plot_level == 'merged_plots':
        # Single figure with multiple subplots
        fig, axes = plt.subplots(1, len(conditions), figsize=(18, 12), sharey=True)
        #legend_handles = set() # List to collect legend handles
        #legend_labels = set() # List to collect legend labels
        legend_handles = [] # List to collect legend handles
        legend_labels = [] # List to collect legend labels
        taxon_colors = {} # Dictionary to store unique colors for each taxon

        
        # Iterate over each condition
        for i, (condition, samples_list) in enumerate(conditions.items()):
            # Select samples corresponding to the condition
            taxo_count = database[[taxo_level] + list(samples_list)]
            taxo_sum = taxo_count.groupby(taxo_level).sum()
            taxo_proportions = taxo_sum.div(taxo_sum.sum(axis=0), axis=1)
            
            # Select top 20 taxonomic categories
            top_20_taxa = taxo_proportions.sum(axis=1).nlargest(20)
            taxo_proportions = taxo_proportions.transpose() # Transposition to have samples in rows and categories in columns
            top_20_taxa = taxo_proportions[top_20_taxa.index]
            other_taxa = 1 - top_20_taxa.sum(axis=1) # Calculation of the sum of the proportions of the remaining categories for each sample
            top_20_taxa['Others'] = other_taxa.values # Addition of the 'Others' line with the proportions of the 'Others' category for each sample
    
            # Plot the barplot for the current condition on the corresponding subplot
            colors = plt.cm.rainbow(range(len(top_20_taxa.columns))) # Generate unique colors for each taxon
            print(colors)
            top_20_taxa.plot(kind='bar', stacked=True, color=colors, ax=axes[i], label=condition, legend=None)
            
            # Get legend handles and labels for each subplot
            handles, labels = axes[i].get_legend_handles_labels()
            legend_handles.extend(handles)
            legend_labels.extend(labels)
            #legend_handles.update(handles)
            #legend_labels.update(labels)
            axes[i].set_ylabel('Abundance')
            axes[i].set_title(condition)
            axes[i].set_xlabel('Samples')
            
        # Store the taxon-color correspondence for this subplot
            for j, taxon in enumerate(top_20_taxa.columns):
                if taxon not in taxon_colors:
                    taxon_colors[taxon] = colors[j]
                    
            top_20_taxa.plot(kind='bar', stacked=True, color=taxon_colors, ax=axes[i], label=condition)

        
        print(taxon_colors)
        # Create a common legend with unique taxon labels and colors
        unique_taxa = list(taxon_colors.keys())
        #unique_colors = [taxon_colors[taxon] for taxon in unique_taxa]
        #for taxon, color in taxon_colors.items():
            #legend_patches.append(mpatches.Patch(color=color, label=taxon))
            
        #plt.legend(legend_handles, unique_taxa, title=taxo_level, loc='best', bbox_to_anchor=(1, 1)) # Create a common legend from the collected handles and labels
        plt.legend(legend_handles, legend_labels, title=taxo_level, loc='best', bbox_to_anchor=(1, 1)) # Create a common legend from the collected handles and labels
        plt.tight_layout()
        format_file = input("Which file format would you like to save the plot ? SVG or PDF ")
        if format_file == 'PDF':
            plt.savefig("barplot_composition_all_samples_condition.pdf", format='pdf', bbox_inches='tight', pad_inches=0.2)
        elif format_file == 'SVG':
            plt.savefig("barplot_composition_all_samples_condition.svg", format='svg', bbox_inches='tight', pad_inches=0.2)
        plt.show()
        exit()

    elif plot_level == 'separate_plots':
        # Barplots for each condition
        for condition, samples_list in conditions.items():
            taxo_count = database[[taxo_level] + list(samples_list)]
            taxo_sum = taxo_count.groupby(taxo_level).sum()
            taxo_proportions = taxo_sum.div(taxo_sum.sum(axis=0), axis=1)

            # Select top 20 taxonomic categories
            top_20_taxa = taxo_proportions.sum(axis=1).nlargest(20)
            taxo_proportions = taxo_proportions.transpose() # Transposition to have samples in rows and categories in columns
            top_20_taxa = taxo_proportions[top_20_taxa.index]
            other_taxa = 1 - top_20_taxa.sum(axis=1) # Calculation of the sum of the proportions of the remaining categories for each sample
            top_20_taxa['Others'] = other_taxa.values # Addition of the 'Others' line with the proportions of the 'Others' category for each sample

            # Plot the barplot for the current condition
            fig, ax = plt.subplots(figsize=(18, 12))
            top_20_taxa.plot(kind='bar', stacked=True, colormap='rainbow', ax=ax, label=condition)
            plt.xlabel('Samples')
            plt.ylabel('Abundance')
            plt.title(f'Abundance of taxonomic category {taxo_level} for condition: {condition}')
            plt.legend(title=taxo_level, loc='best', bbox_to_anchor=(1, 1))
            plt.xticks(rotation=90, ha='right')
            plt.tight_layout()
            format_file = input("Which file format would you like to save the plot ? SVG or PDF ")
            if format_file == 'PDF':
                plt.savefig(f"barplot_composition_{condition}.pdf", format='pdf', bbox_inches='tight', pad_inches=0.2)
            elif format_file == 'SVG':
                plt.savefig(f"barplot_composition_{condition}.svg", format='svg', bbox_inches='tight', pad_inches=0.2)
            plt.show()
    else:
        print("Plots representation not supported.")
        
        

def barplot_all_samples_two_conditions(database, condition, asv_info):
    # Extract column names of samples
    sequence_index = list(database.columns).index("sequence")
    samples_columns = database.columns[sequence_index + 1:]
    
    taxo_level = input("What taxonomic level do you want ? \nFor assignment with PR2 : Domain;Supergroup;Division;Subdivision;Class;Order;Family;Genus;Species \nFor assignment with Silva : Kingdom;Phylum;Class;Order;Family;Genus \nYour choice : ")
    
    conditions=two_conditions_struct(asv_info, condition) # Get the list of samples associated with the condition
    
    plot_level = input("Do you want a merged chart with all conditions or separate charts for each condition ? merged_plots / separate_plots ")
    
    if plot_level == 'merged_plots':
        # Initialize a dictionary to store summed abundances for each condition
        condition_abundances = {cond: 0 for cond in conditions}
        
        # Iterate over each conditions
        for i, ((cond1, cond2), samples_list) in enumerate(conditions.items()):
            # Select samples corresponding to the condition
            taxo_count = database[[taxo_level] + list(samples_list)]
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

        color_map = sns.color_palette("tab20", len(all_taxons))
        
        for i, (condition_name, taxon_abundances) in enumerate(condition_abundances.items()):
            top_20_taxa = taxon_abundances.sort_values(ascending=False)[:20]
            other_taxa = 1 - top_20_taxa.sum(axis=0)
            top_20_taxa.loc['Others'] = other_taxa # Addition of the 'Others' line with the proportions of the 'Others' category for each condition
            
            # Initialize y_offset to keep track of the vertical position of the next bar
            y_offset = 0
            for j, (taxon, abundance) in enumerate(top_20_taxa.items()):
                if taxon in top_20_taxa:
                    if taxon != 'Others':
                        color = color_map[j] # Assign color based on taxon
                    else:
                        color = 'yellow' # Assign color for 'Others' categorie
                    plt.bar((f"{condition_name[0]} - {condition_name[1]}"), abundance, color=color, bottom=y_offset)  # Plot each taxon with assigned color
                    y_offset += abundance
                    if i == 0:  # Only add handles from the first condition to avoid duplicates
                        legend_handles.append(taxon)
            
            
            plt.text(i, -0.05, condition_name[0], ha='center', va='bottom', fontsize=12, fontweight='bold')
            plt.text(i, 1.04, condition_name[1], ha='center', va='top', fontsize=12, fontweight='bold')

                    
        plt.xlabel('Condition')
        plt.ylabel('Total Abundance')
        plt.title('Total Abundance by Condition')
        plt.xticks(rotation=90)
        plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        plt.legend(legend_handles, title=taxo_level, loc='best', bbox_to_anchor=(1, 1))
        plt.tight_layout()
        
        format_file = input("Which file format would you like to save the plot ? SVG or PDF ")
        if format_file == 'PDF':
            plt.savefig("barplot_total_abundance_by_two_condition.pdf", format='pdf', bbox_inches='tight', pad_inches=0.2)
        elif format_file == 'SVG':
            plt.savefig("barplot_total_abundance_by_two_condition.svg", format='svg', bbox_inches='tight', pad_inches=0.2)
        plt.show()

    elif plot_level == 'separate_plots':
        # Barplots for each condition
        for condition, samples_list in conditions.items():
            taxo_count = database[[taxo_level] + list(samples_list)]
            taxo_sum = taxo_count.groupby(taxo_level).sum()
            taxo_proportions = taxo_sum.div(taxo_sum.sum(axis=0), axis=1)

            # Select top 20 taxonomic categories
            top_20_taxa = taxo_proportions.sum(axis=1).nlargest(20)
            taxo_proportions = taxo_proportions.transpose() # Transposition to have samples in rows and categories in columns
            top_20_taxa = taxo_proportions[top_20_taxa.index]
            other_taxa = 1 - top_20_taxa.sum(axis=1) # Calculation of the sum of the proportions of the remaining categories for each sample
            top_20_taxa['Others'] = other_taxa.values # Addition of the 'Others' line with the proportions of the 'Others' category for each sample

            # Plot the barplot for the current condition
            fig, ax = plt.subplots(figsize=(18, 12))
            top_20_taxa.plot(kind='bar', stacked=True, colormap='rainbow', ax=ax, label=condition)
            plt.xlabel('Samples')
            plt.ylabel('Abundance')
            plt.title(f'Abundance of taxonomic category {taxo_level} for condition: {condition}')
            plt.legend(title=taxo_level, loc='best', bbox_to_anchor=(1, 1))
            plt.xticks(rotation=90, ha='right')
            plt.tight_layout()
            format_file = input("Which file format would you like to save the plot ? SVG or PDF ")
            if format_file == 'PDF':
                plt.savefig(f"barplot_composition_{condition}.pdf", format='pdf', bbox_inches='tight', pad_inches=0.2)
            elif format_file == 'SVG':
                plt.savefig(f"barplot_composition_{condition}.svg", format='svg', bbox_inches='tight', pad_inches=0.2)
            plt.show()
    else:
        print("Plots representation not supported.")
        
        
        

def heatmap(database): # Function to create a heatmap of taxonomic composition for all samples
    # Extract column names of samples
    sequence_index = list(database.columns).index("sequence")
    samples_columns = database.columns[sequence_index + 1:]
    
    taxo_level = input("What taxonomic level do you want ? \nFor assignment with PR2 : Domain;Supergroup;Division;Subdivision;Class;Order;Family;Genus;Species \nFor assignment with Silva : Kingdom;Phylum;Class;Order;Family;Genus \nYour choice : ")
    
    # Calculate proportions of each taxonomic category for each sample
    taxo_count = database[[taxo_level] + list(samples_columns)]
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
    format_file = input("Which file format would you like to save the plot ? SVG or PDF ")
    if format_file == 'PDF':
        plt.savefig("heatmap_proportion.pdf", format='pdf', bbox_inches='tight', pad_inches=0.2)
    elif format_file == 'SVG':
        plt.savefig("heatmap_proportion.svg", format='svg', bbox_inches='tight', pad_inches=0.2)
    plt.show()
    exit()


def piechart(database): # Function to create a pie chart of taxonomic composition for all samples
    # Extract column names of samples
    sequence_index = list(database.columns).index("sequence")
    samples_columns = database.columns[sequence_index + 1:]
    
    taxo_level = input("What taxonomic level do you want ? \nFor assignment with PR2 : Domain;Supergroup;Division;Subdivision;Class;Order;Family;Genus;Species \nFor assignment with Silva : Kingdom;Phylum;Class;Order;Family;Genus \n ASVNumber \nYour choice : ")
    
    taxo_count = database[[taxo_level] + list(samples_columns)] # Select columns related to the chosen taxonomic level
    taxo_sum = taxo_count.groupby(taxo_level).sum() # Group by the chosen taxonomic level and sum the counts
    total_counts = taxo_sum.sum(axis=1) # Sum the counts across all samples
    sorted_total_counts = total_counts.sort_values(ascending=False)[:10]  # Sort the top 10 total counts in descending order
    other_taxa = total_counts.sort_values(ascending=False)[10:] # Sort the rest
    
    # Calculate the values for the pie chart
    values = sorted_total_counts.values
    other_percent = other_taxa.sum()
    other_series = pd.Series(other_percent, index=['Others'])
    filtered_families = pd.Series(values, index=sorted_total_counts.index)
    filtered_families = pd.concat([filtered_families, other_series]) # Merge dataframe of top 20 and 'Others' categories

    # Pie chart of taxonomic composition
    plt.figure(figsize=(12, 8))
    plt.pie(filtered_families, labels=filtered_families.index, autopct='%1.1f%%', startangle=140)
    plt.title(f'Proportion of taxonomic category {taxo_level} in all samples')
    plt.axis('equal')
    format_file = input("Which file format would you like to save the plot ? SVG or PDF ")
    if format_file == 'PDF':
        plt.savefig("piechart.pdf", format='pdf', pad_inches=0.2)
    elif format_file == 'SVG':
        plt.savefig("piechart.svg", format='svg', pad_inches=0.2)
    plt.show()
    exit()
        
    
def venn_diagram(database): # Function to create a Venn diagram of taxonomic overlap between samples
    # Extract column names of samples
    sequence_index = list(database.columns).index("sequence")
    samples_columns = database.columns[sequence_index + 1:] # Sélection des colonnes correspondantes aux échantillons (après la colonne "sequence")
    
    taxo_level = input("What taxonomic level do you want ? \nFor assignment with PR2 : Domain;Supergroup;Division;Subdivision;Class;Order;Family;Genus;Species \nFor assignment with Silva : Kingdom;Phylum;Class;Order;Family;Genus \nNumber of ASV; \nYour choice : ")
    if taxo_level == 'ASV':
        taxo_level="ASVNumber"
        
    # Loop through each sample column to identify taxa present in each sample
    taxo_samples = [set(database[database[colonne].notna() & (database[colonne] > 0)][taxo_level]) for colonne in samples_columns]
    
    # Choose the appropriate Venn diagram based on the number of samples
    if len(samples_columns) == 2:
        venn2(taxo_samples, (samples_columns[0], samples_columns[1]))
    elif len(samples_columns) == 3:
        venn3(taxo_samples, (samples_columns[0], samples_columns[1], samples_columns[2]))
    #elif len(samples_columns) > 3:
        #upset
    else:
        print("Number of samples not currently supported.")
    
    if taxo_level == 'ASV':
        plt.title(f'Venn diagram of ASV in samples')
    else:
        plt.title(f'Venn diagram of taxonomic category {taxo_level} in samples')
    plt.tight_layout()
    format_file = input("Which file format would you like to save the plot ? SVG or PDF ")
    if format_file == 'PDF':
        plt.savefig("venn_diagram.pdf", format='pdf', pad_inches=0.5)
    elif format_file == 'SVG':
        plt.savefig("venn_diagram.svg", format='svg', pad_inches=0.5)
    plt.show()
    exit()
        
        
def composition_diagram(database):  # Function to create a composition diagram showing the distribution of taxonomic categories in samples
    # Extract column names of samples
    sequence_index = list(database.columns).index("sequence")
    samples_columns = database.columns[sequence_index + 1:]

    taxo_level = input("What taxonomic level do you want ? \nFor assignment with PR2 : Domain;Supergroup;Division;Subdivision;Class;Order;Family;Genus;Species \nFor assignment with Silva : Kingdom;Phylum;Class;Order;Family;Genus\n \nYour choice : ")
            
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
    format_file = input("Which file format would you like to save the plot ? SVG or PDF ")
    if format_file == 'PDF':
        plt.savefig("diagram_taxo_composition.pdf", format='pdf', bbox_inches='tight', pad_inches=0.2)
    elif format_file == 'SVG':
        plt.savefig("diagram_taxo_composition.svg", format='svg', bbox_inches='tight', pad_inches=0.2)
    plt.show()
    exit()
    
    
def network(asv_taxo): # Function to create a network graph based on taxonomic similarity
        taxo_level = input("What taxonomic level do you want ? \nFor assignment with PR2 : Domain;Supergroup;Division;Subdivision;Class;Order;Family;Genus;Species \nFor assignment with Silva : Kingdom;Phylum;Class;Order;Family;Genus \nYour choice : ")
        G = nx.Graph() # Undirected graph
       
        # Add nodes corresponding to ASV numbers
        for asv in asv_taxo["ASVNumber"]:
            G.add_node(asv)

        # Add edges between ASVs based on similarity for the specified taxonomic level
        for i, row1 in asv_taxo.iterrows():
            for j, row2 in asv_taxo.iterrows():
                if i < j:
                    if row1[taxo_level] == row2[taxo_level] and row1["ASVNumber"] != row2["ASVNumber"]:
                        G.add_edge(row1["ASVNumber"], row2["ASVNumber"])

        nx.draw(G, with_labels=True)
        format_file = input("Which file format would you like to save the plot ? SVG or PDF ")
        if format_file == 'PDF':
            plt.savefig("network.pdf", format='pdf', pad_inches=0.2)
        elif format_file == 'SVG':
            plt.savefig("network.svg", format='svg', pad_inches=0.2)
        plt.show()
        exit()
            


####################
####### MAIN #######
####################

if __name__ == "__main__":
    asv_info=sys.argv[1] # Path to "asv.csv" file, output from dada2 [mandatory]
    asv_taxo=sys.argv[2] # Path to "taxo.csv" file, output from dada2 [mandatory]
    database=sys.argv[3] # Path to "database.csv" file, output from dada2 [mandatory]
    
    if len(sys.argv) >= 5:
        condition = sys.argv[4] # Path to the file (csv or excel) containing the names of the samples and the conditions to which the samples belong
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
    scaler_choice = input("Which scaler would you like for data normalization ? MinMax / StandardScaler : ")
    if scaler_choice == 'MinMax':
        scaler = MinMaxScaler()
    elif scaler_choice == 'StandardScaler':
        scaler = StandardScaler()
    else :
        print("Scaler choice not supported.")
     # Normalization of asv_info
    asv_info_normalized = asv_info.copy()
    columns_numeric = asv_info.columns[1:]
    asv_info_normalized[columns_numeric] = scaler.fit_transform(asv_info_normalized[columns_numeric])
    # Normalization of database
    database_normalized = database.copy()
    sequence_index = list(database.columns).index("sequence")
    columns_numeric = database.columns[sequence_index + 1:]
    database_normalized[columns_numeric] = scaler.fit_transform(database_normalized[columns_numeric])
    
    
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
                statistical_test_alpha(asv_info_normalized, condition)
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
                df_counts, beta_diversity = beta_diversity_function(asv_info_normalized)
                beta_diversity_all(beta_diversity)
            elif analyse_beta == 'statistical_test':
                df_counts, beta_diversity = beta_diversity_function(asv_info_normalized)
                statistical_test_beta(df_counts, beta_diversity)
            elif analyse_beta == 'beta_diversity_graph':
                df_counts, beta_diversity = beta_diversity_function(asv_info_normalized)
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
                print("Type of representation not supported.")
                exit()

        else:
            print("Analysis not supported.")
            
        keep_going = input("Do you want to run another analysis ? (yes / no) ")

