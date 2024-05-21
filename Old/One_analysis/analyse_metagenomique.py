import sys
import pandas as pd
import skbio
from skbio.diversity import beta_diversity
from skbio.stats.ordination import pcoa
from skbio.stats.distance import permanova
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
import seaborn as sns
import numpy as np
from scipy.stats import f_oneway
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import ast
import networkx as nx
from upsetplot import UpSet
from upsetplot import plot
from upsetplot import generate_counts
from itertools import combinations


#####################
## Alpha Diversity ##
#####################

def alpha_diversity_one(asv_info):
    print("Available samples :")
    for column in asv_info.columns[2:]:
        print(column)
    sample_alpha = input("Which sample do you want alpha diversity ? : ")
    counts = asv_info[sample_alpha]
    alpha_index = input("Which alpha diversity index do you want to calculate ? (shannon / simpson): ")
    if alpha_index == 'shannon':
        alpha_diversity = skbio.diversity.alpha.shannon(counts)
        print("Alpha diversity : ", alpha_diversity)
    elif alpha_index == 'simpson':
        alpha_diversity = skbio.diversity.alpha.simpson(counts)
        print("Alpha diversity : ", alpha_diversity)
    else:
        print("Alpha diversity index not supported.")
        exit()
    
   
def alpha_diversity_all(asv_info):
    alpha_diversity_all = {}
    alpha_index = input("Which alpha diversity index do you want to calculate ? (shannon / simpson): ")
    if alpha_index == 'shannon':
        for column in asv_info.columns[2:]:  # Ignorer les deux premières colonnes (ASVNumber et sequence)
            counts = asv_info[column]
            alpha_diversity_all[column] = {}
            alpha_diversity_all[column]['shannon'] = skbio.diversity.alpha.shannon(counts)
        # Affichage des diversités alpha pour tous les échantillons
        for column, diversity in alpha_diversity_all.items():
            print("-- Alpha diversity for sample ", column, " : ")
            print("Shannon:", diversity['shannon'])
    elif alpha_index == 'simpson':
        for column in asv_info.columns[2:]:
            counts = asv_info[column]
            alpha_diversity_all[column] = {}
            alpha_diversity_all[column]['simpson'] = skbio.diversity.alpha.simpson(counts)
        # Affichage des diversités alpha pour tous les échantillons
        for column, diversity in alpha_diversity_all.items():
            print("-- Alpha diversity for sample ", column, " : ")
            print("Simpson:", diversity['simpson'])
    else:
        print("Alpha diversity index not supported.")
        exit()
    
    
def statistical_test_alpha(asv_info):
    type_test = input("What type of statistical test would you like to perform? Comparison of each variable (anova) or of two variables together (post_hoc_tuckey) : ")
    if type_test == 'anova':
        counts = [asv_info[col] for col in asv_info.columns[2:]]  # Ignorer les deux premières colonnes (ASVNumber et sequence)
        f_statistic, p_value = f_oneway(*counts)
        print("ANOVA test result : ")
        print("Statistic F : ", f_statistic)
        print("p-value : ", p_value)
    elif type_test == 'post_hoc_tuckey': #95% test
        flattened_data = [asv_info[col].dropna() for col in asv_info.columns[2:]]
        group_labels = asv_info.columns[2:]
        for data in flattened_data:
            tukey_result = pairwise_tukeyhsd(np.concatenate(flattened_data), np.repeat(group_labels, len(data)))
        print("Result of Tukey test : ")
        print(tukey_result)
    else:
        print("Type of test not supported.")
        exit()
    


def alpha_graph(asv_info):
    alpha_index = input("Which alpha diversity index would you like ? (shannon / simpson / both) : ")
    if alpha_index == 'both':
        alpha_index = ['shannon', 'simpson']
    elif alpha_index not in ['shannon', 'simpson']:
        print("Alpha diversity index not supported.")
        exit()
    else:
        alpha_index = [alpha_index]
        
    alpha_diversity = {}
    for column in asv_info.columns[2:]:  # Ignorer les deux premières colonnes (ASVNumber et sequence)
        counts = asv_info[column]
        for alpha_index in alpha_index:
            if alpha_index == 'shannon':
                alpha_diversity.setdefault('shannon', {})[column] = skbio.diversity.alpha.shannon(counts)
            elif alpha_index == 'simpson':
                alpha_diversity.setdefault('simpson', {})[column] = skbio.diversity.alpha.simpson(counts)

    # Création de la figure avec deux sous-graphiques
    fig, axs = plt.subplots(1, len(alpha_index), figsize=(14, 6))
    if len(alpha_index) == 1:
        axs = [axs]  # Si un seul indice est choisi, axs est une liste simple au lieu d'une liste de listes
    
    for i, alpha_index in enumerate(alpha_index):
        # Tracer le graphique pour chaque indice de diversité alpha
        axs[i].scatter(alpha_diversity[alpha_index].keys(), alpha_diversity[alpha_index].values(), color='blue', marker='o')
        axs[i].set_xlabel('Samples')
        axs[i].set_ylabel(f'{alpha_index.capitalize()} alpha diversity')
        axs[i].set_title(f'{alpha_index.capitalize()} alpha diversity according to samples')
        axs[i].tick_params(axis='x', rotation=45, labelsize=8)
    plt.tight_layout()
    plt.savefig("alpha_diversity.pdf", format='pdf')
    plt.show()
    


####################
## Beta Diversity ##
####################

def beta_diversity_all(asv_info):
    count_table = {}
    for column in asv_info.columns[2:]:  # Ignorer les deux premières colonnes (ASVNumber et sequence)
        asv_counts = {}
        for index, row in asv_info.iterrows():
            count = row[column] # Récupérer le comptage de cet ASV dans cet échantillon
            asv_counts[index] = count
        count_table[column] = asv_counts # Ajouter les comptages de cet échantillon à la table
    df_counts = pd.DataFrame(count_table).T.fillna(0)  # Transposer pour que les échantillons soient les lignes
        
    print(count_table)
    print("df_counts")
    print(df_counts)
        
    beta_index = input("Which alpha diversity indices would you like ? : (euclidean/cityblock/braycurtis/canberra/chebyshev/correlation/cosine/dice/hamming/jaccard/kulsinski/mahalanobis/manhattan/matching/minkowski/rogerstanimoto/russellrao/seuclidean/sokalmichener/sokalsneath/sqeuclidean/yule)")
    if beta_index not in ['braycurtis', 'canberra', 'chebyshev', 'cityblock', 'correlation', 'cosine', 'dice', 'euclidean', 'hamming', 'jaccard', 'kulsinski', 'mahalanobis', 'manhattan', 'matching', 'minkowski', 'rogerstanimoto', 'russellrao', 'seuclidean', 'sokalmichener', 'sokalsneath', 'sqeuclidean', 'yule']:
        print("Beta diversity index not supported. The available index are : braycurtis, canberra, chebyshev, cityblock, correlation, cosine, dice, euclidean, hamming, jaccard, kulsinski, mahalanobis, manhattan, matching, minkowski, rogerstanimoto, russellrao, seuclidean, sokalmichener, sokalsneath, sqeuclidean, yule ")
        exit()
    beta_diversity = skbio.diversity.beta_diversity(beta_index, df_counts)
    print("-- Beta diversity :")
    print(beta_diversity)
    
    return df_counts, beta_diversity
        
def statistical_test_beta(df_counts, beta_diversity):
    groups = [nom.split('_')[0] for nom in df_counts.index]   # Extraire les groupes à partir des noms des échantillons
    print(groups)
    permanova_result = permanova(distance_matrix=beta_diversity, grouping=groups)
    print("-- Permanova test results")
    print("Permanova statistic : ", permanova_result['test statistic'])
    print("p-value : ", permanova_result['p-value'])
    print("Number of permutations : ", permanova_result['number of permutations'])
    
    
def beta_diversity_graph(beta_diversity):
    beta_representation = input("Which beta diversity representation would you like ? (heatmap / PCoA): ")
    if beta_representation == 'heatmap':
        #beta_diversity=pd.DataFrame(beta_diversity)
        plt.figure(figsize=(10, 8))
        sns.heatmap(beta_diversity.to_data_frame(), cmap="viridis", annot=True, fmt=".2f", linewidths=.5)
        plt.title("Beta diversity heatmap")
        plt.xlabel("Samples")
        plt.ylabel("Samples")
        plt.savefig("heatmap_beta_diversity.pdf", format='pdf')
        plt.show()
    elif beta_representation == 'PCoA':
        print(beta_diversity)
        pcoa_results = pcoa(beta_diversity)
        if pcoa_results.samples.shape[1] >= 3:
            beta_dimension = input("Which dimension would you like? (3D / 2D): ")
            if beta_dimension == '3D':
            # 3 dimentions
                pcoa_results.plot(cmap='viridis', title='PCoA Plot')
                plt.savefig("PCoA_3D.pdf", format='pdf')
                plt.show()
            elif beta_dimension == '2D':
            # 2 dimensions, coloré par échantillons
                plt.scatter(pcoa_results.samples.iloc[:, 0], pcoa_results.samples.iloc[:, 1], c=range(len(pcoa_results.samples)), cmap='viridis')
                plt.title('PCoA Plot')
                plt.xlabel('PCoA 1')
                plt.ylabel('PCoA 2')
                plt.colorbar(label='Samples')
                plt.savefig("PCoA_2D.pdf", format='pdf')
                plt.show()
            else:
                print("Dimension not supported.")
        else:
            print("There are not enough dimensions available to trace the PCoA.")
    else:
        print("Beta diversity representation not supported.")
        exit()
    


#################################
## Graph taxonomic composition ##
#################################


def barplot_one_sample(database):
     #barplot de la composition du niveau taxonomique dans l'échantillon choisi
    print("Available samples :")
    available_samples = list(database.columns[database.columns.get_loc('sequence')+1:])
    for sample in available_samples:
        print(sample)
    sample = input("Which sample would you like to analyze ? ")
    taxo_level = input("What taxonomic level do you want ? \nFor assignment with PR2 : Domain;Supergroup;Division;Subdivision;Class;Order;Family;Genus;Species \nFor assignment with Silva : Kingdom;Phylum;Class;Order;Family;Genus \nYour choice : ")

    taxo_count = database[[taxo_level] + [sample]]
    taxo_sum = taxo_count.groupby(taxo_level).sum()
    taxo_sum[sample] = taxo_sum[sample].astype(float)
    taxo_proportions = taxo_sum / taxo_sum.sum()

    # Fusionner les taxons avec des proportions inférieures à 0.001 dans une catégorie "Autres"
    taxo_others = taxo_proportions[taxo_proportions[sample] < 0.001]
    proportion_others = taxo_others.sum()
    taxo_proportions_without_others = taxo_proportions[taxo_proportions[sample] >= 0.001]
    taxo_proportions_finals = pd.concat([taxo_proportions_without_others, pd.DataFrame([proportion_others], index=['Autres'], columns=[sample])])

    plt.figure(figsize=(10, 6))
    taxo_proportions_finals[sample].sort_values(ascending=False).plot(kind='bar')
    plt.xlabel(taxo_level)
    plt.ylabel('Proportion')
    plt.title(f'Proportion of {taxo_level} taxonomic category in sample {sample}')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig("barplot_one_sample.pdf", format='pdf')
    plt.show()

    
def barplot_all_samples(database): #barplot de la composition du niveau taxonomique dans les différents échantillons
    # Trouver l'indice de la colonne "sequence"
    sequence_index = list(database.columns).index("sequence")
    # Sélectionner les colonnes correspondant aux échantillons (après la colonne "sequence")
    samples_columns = database.columns[sequence_index + 1:]
    taxo_level = input("What taxonomic level do you want ? \nFor assignment with PR2 : Domain;Supergroup;Division;Subdivision;Class;Order;Family;Genus;Species \nFor assignment with Silva : Kingdom;Phylum;Class;Order;Family;Genus \nYour choice : ")
    taxo_count = database[[taxo_level] + list(samples_columns)]
    taxo_sum = taxo_count.groupby(taxo_level).sum()
    taxo_proportions = taxo_sum.div(taxo_sum.sum(axis=0), axis=1)
    taxo_proportions_transpose = taxo_proportions.transpose()
        
    plt.figure(figsize=(40, 40))
    taxo_proportions_transpose.plot(kind='bar', stacked=True, colormap='rainbow')
    plt.xlabel('Samples')
    plt.ylabel('Abundance')
    plt.title(f'Abundance of taxonomic category {taxo_level} for all samples')
    plt.xticks(rotation=45, ha='right')
    plt.savefig("barplot_composition_all_samples.pdf", format='pdf')
    plt.tight_layout()
    plt.show()
        
def heatmap(database): # Heatmap de la composition du niveau taxonomique dans les différents échantillons
    sequence_index = list(database.columns).index("sequence")
    samples_columns = database.columns[sequence_index + 1:]
    taxo_level = input("What taxonomic level do you want ? \nFor assignment with PR2 : Domain;Supergroup;Division;Subdivision;Class;Order;Family;Genus;Species \nFor assignment with Silva : Kingdom;Phylum;Class;Order;Family;Genus \nYour choice : ")
    taxo_count = database[[taxo_level] + list(samples_columns)]
    taxo_sum = taxo_count.groupby(taxo_level).sum()
    taxo_proportions = taxo_sum.div(taxo_sum.sum(axis=0), axis=1)
    taxo_proportions_heatmap = taxo_proportions.transpose()

    plt.figure(figsize=(12, 8))
    sns.heatmap(taxo_proportions_heatmap, cmap='rainbow')
    plt.xlabel(taxo_level)
    plt.ylabel('Samples')
    plt.title(f'Heatmap of the proportion of taxonomic category {taxo_level} in each sample')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig("heatmap_proportion.pdf", format='pdf')
    plt.show()


def piechart(database): # Diagramme circulaire de la composition du niveau taxonomique choisi dans tous les échantillons confondus
    sequence_index = list(database.columns).index("sequence")
    samples_columns = database.columns[sequence_index + 1:]
    taxo_level = input("What taxonomic level do you want ? \nFor assignment with PR2 : Domain;Supergroup;Division;Subdivision;Class;Order;Family;Genus;Species \nFor assignment with Silva : Kingdom;Phylum;Class;Order;Family;Genus \nYour choice : ")
    taxo_count = database[[taxo_level] + list(samples_columns)]
    taxo_sum = taxo_count.groupby(taxo_level).sum()
    total_counts = taxo_sum.sum(axis=1)
    sorted_total_counts = total_counts.sort_values(ascending=False)
    values = sorted_total_counts.values
    values_sum = values.sum()
        
    # Fusionner les taxons avec des proportions inférieures à 0.01 dans une catégorie "Autres" afin de rendre le graphique plus lisible
    percent = 0.01*values_sum
    filtered_families = values[values >= percent]
    other_percent = values[values < percent].sum()
    other_series = pd.Series(other_percent, index=['Autres'])
    filtered_families = pd.Series(filtered_families, index=sorted_total_counts.index[sorted_total_counts >= percent])
    filtered_families = pd.concat([filtered_families, other_series])
    print(filtered_families)

    plt.figure(figsize=(12, 8))
    plt.pie(filtered_families, labels=filtered_families.index, autopct='%1.1f%%', startangle=140) #, rotatelabels=True
    plt.title(f'Proportion de la catégorie taxonomique {taxo_level} dans tous les échantillons')
    plt.axis('equal')
    plt.savefig("piechart.pdf", format='pdf')
    plt.show()
        
    
def venn_diagram(database):
    sequence_index = list(database.columns).index("sequence")
    samples_columns = database.columns[sequence_index + 1:] # Sélection des colonnes correspondantes aux échantillons (après la colonne "sequence")
    taxo_level = input("What taxonomic level do you want ? \nFor assignment with PR2 : Domain;Supergroup;Division;Subdivision;Class;Order;Family;Genus;Species \nFor assignment with Silva : Kingdom;Phylum;Class;Order;Family;Genus \nYour choice : ")
    for col in samples_columns:
        taxo_samples = set(database[database[col].notna() & (database[col] > 0)][taxo_level]) # Identifier les taxons présents dans chaque échantillon
    if len(samples_columns) == 2:
        venn2(taxo_samples, (samples_columns[0], samples_columns[1]))
    elif len(samples_columns) == 3:
        venn3(taxo_samples, (samples_columns[0], samples_columns[1], samples_columns[2]))
    #elif len(samples_columns) > 3:
        #upset
        
    else:
        print("Number of samples not currently supported.")
        
    plt.title(f'Venn diagram of taxonomic category {taxo_level} in samples')
    plt.savefig("venn.pdf", format='pdf')
    plt.show()
        
        
def composition_diagram(database):
    sequence_index = list(database.columns).index("sequence")
    samples_columns = database.columns[sequence_index + 1:] # Sélection des colonnes correspondantes aux échantillons (après la colonne "sequence")
    taxo_level = input("What taxonomic level do you want ? \nFor assignment with PR2 : Domain;Supergroup;Division;Subdivision;Class;Order;Family;Genus;Species \nFor assignment with Silva : Kingdom;Phylum;Class;Order;Family;Genus \nYour choice : ")
            
    # Identifier les taxons présentes dans chaque échantillon et sa frequence d'apparition
    for col in samples_columns:
        taxo_samples = set(database[database[col].notna() & (database[col] > 0)][taxo_level])
    print(taxo_samples)
    taxo_samples = {}
    for col in samples_columns:
        taxo_samples[col] = database[database[col].notna() & (database[col] > 0)][taxo_level].value_counts(normalize=True)
    data = pd.DataFrame(taxo_samples)

    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(111)
    data.plot(kind='bar', ax=ax)
    plt.xlabel(taxo_level)
    plt.ylabel('Frequency in samples')
    plt.title(f'Diagram of {taxo_level} taxonomic category distribution in samples')
    plt.legend(title="Samples", loc="upper right", bbox_to_anchor=(1.2, 1))
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig("diagram_taxo_composition.pdf", format='pdf')
    plt.show()
    
    
def network(asv_taxo):
        taxo_level = input("What taxonomic level do you want ? \nFor assignment with PR2 : Domain;Supergroup;Division;Subdivision;Class;Order;Family;Genus;Species \nFor assignment with Silva : Kingdom;Phylum;Class;Order;Family;Genus \nYour choice : ")
        G = nx.Graph() # Graphe non orienté
        for asv in asv_taxo["ASVNumber"]:
            G.add_node(asv) # Ajouter les nœuds correspondant aux numéros ASV

        # Ajouter les arêtes entre les ASV en fonction de la similarité pour le taxon spécifié
        for i, row1 in asv_taxo.iterrows():
            for j, row2 in asv_taxo.iterrows():
                if i < j:
                    if row1[taxo_level] == row2[taxo_level] and row1["ASVNumber"] != row2["ASVNumber"]:
                        G.add_edge(row1["ASVNumber"], row2["ASVNumber"])

        nx.draw(G, with_labels=True)
        plt.show()
            



####################
####### MAIN #######
####################

if __name__ == "__main__":
    asv_info=sys.argv[1] # Path to "asv.csv" file, output from dada2
    asv_taxo=sys.argv[2] # Path to "taxo.csv" file, output from dada2
    database=sys.argv[3] # Path to "database.csv" file, output from dada2
    asv_info = pd.read_csv(asv_info)
    asv_taxo = pd.read_csv(asv_taxo)
    database = pd.read_csv(database)
    
    keep_going = 'yes'
    while keep_going.lower() == 'yes':
        # Demander à l'utilisateur ce qu'il veut étudier
        analyse = input("What would you like to study? (alpha_diversity / beta_diversity / graph_taxo_composition):")
        
        if analyse == 'alpha_diversity':
            analyse_alpha = input("What you want to know about alpha diversity ? (diversity_one_sample / diversity_all_samples / statistical_test / alpha_diversity_graph) : ")
            if analyse_alpha == 'diversity_one_sample':
                alpha_diversity_one(asv_info)
            elif analyse_alpha == 'diversity_all_samples':
                alpha_diversity_all(asv_info)
            elif analyse_alpha == 'statistical_test':
                statistical_test_alpha(asv_info)
            elif analyse_alpha == 'alpha_diversity_graph':
                alpha_graph(asv_info)
            else:
                print("Analysis not supported.")
                exit()
            
        elif analyse == 'beta_diversity':
            analyse_beta = input("What you want to know about beta diversity ? (beta_diversity / statistical_test / beta_diversity_graph) : ")
            if analyse_beta == 'beta_diversity':
                beta_diversity_all(asv_info)
            elif analyse_beta == 'statistical_test':
                df_counts, beta_diversity = beta_diversity_all(asv_info)
                statistical_test_beta(df_counts, beta_diversity)
            elif analyse_beta == 'beta_diversity_graph':
                df_counts, beta_diversity = beta_diversity_all(asv_info)
                beta_diversity_graph(beta_diversity)
            else:
                print("Analysis not supported.")
                exit()
            
        elif analyse == 'graph_taxo_composition':
            composition_level = input("What kind of representation would you like ? (barplot_one_sample / barplot_all_samples / heatmap / piechart / venn_diagram / diagram_taxo_composition /network ): ")
            if composition_level == 'barplot_one_sample':
                barplot_one_sample(database)
            if composition_level == 'barplot_all_samples':
                barplot_all_samples(database)
            if composition_level == 'heatmap':
                heatmap(database)
            if composition_level == 'piechart':
                piechart(database)
            if composition_level == 'venn_diagram':
                venn_diagram(database)
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

