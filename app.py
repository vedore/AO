import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt

from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.model_selection import train_test_split


def gene_assignment(gene_df):
    new_gene_def = []

    for gene in gene_df:
        gene_disp = []
        lines = gene.split("///")

        for line in lines:
            new_line = line.strip().split("//")

            if len(new_line) == 5:
                gene_disp.append(new_line[1].strip())

        new_gene_def.append(gene_disp)

    return new_gene_def


def join_column_dispersion(dt, file_dispersion):
    df_file_dispersion = pd.read_csv(f'./data_files/GSMFiles/{file_dispersion}', header=None, delimiter="\t")
    df_file_dispersion.columns = ['id', 'disp']

    dt.loc[:, 'disp'] = pd.Series([-1.0 for x in range(len(dt.index))])

    for x in range(0, len(df_file_dispersion)):
        id_value = df_file_dispersion['id'][x]
        disp_value = df_file_dispersion['disp'][x]

        dt.loc[dt.id == id_value, 'disp'] = disp_value

    return dt


def scale01(x):
    return (x - x.min()) / (x.max() - x.min())


def clustering_ips_activin_esc_activin():
    df_ips_activin_26 = pd.read_csv('./data_files/GSMFiles/GSM958226-tbl-1.txt', header=None, delimiter="\t")
    df_ips_activin_27 = pd.read_csv('./data_files/GSMFiles/GSM958227-tbl-1.txt', header=None, delimiter="\t")
    df_esc_activin_20 = pd.read_csv('./data_files/GSMFiles/GSM958220-tbl-1.txt', header=None, delimiter="\t")
    df_esc_activin_21 = pd.read_csv('./data_files/GSMFiles/GSM958221-tbl-1.txt', header=None, delimiter="\t")

    df_ips_activin_26.columns = ['id', 'disp']
    df_ips_activin_27.columns = ['id', 'disp']
    df_esc_activin_20.columns = ['id', 'disp']
    df_esc_activin_21.columns = ['id', 'disp']

    df = pd.DataFrame()
    df['id'] = df_ips_activin_26['id']
    df['ips_26'] = df_ips_activin_26['disp']
    df['ips_27'] = df_ips_activin_27['disp']
    df['esc_20'] = df_esc_activin_20['disp']
    df['esc_21'] = df_esc_activin_21['disp']

    # for var in range(0, d):
    #     data[:, var] = scale01(data[:, var])

    selected_columns = ['ips_26', 'ips_27', 'esc_20', 'esc_21']
    data = df[selected_columns]

    kmeans = KMeans(n_clusters=3, random_state=0, n_init='auto')
    kmeans.fit(data)

    cluster_labels = kmeans.labels_

    pca = PCA(n_components=2)
    reduced_data = pca.fit_transform(data)

    print("ok")

    plt.scatter(reduced_data[:, 0], reduced_data[:, 1], c=cluster_labels, cmap='viridis', alpha=0.5, label='Data Points')
    plt.title('KMeans Clustering (PCA Reduced)')
    plt.xlabel('Principal Component 1')
    plt.ylabel('Principal Component 2')
    plt.legend()
    plt.show()


clustering_ips_activin_esc_activin()

# df_all_values = pd.read_csv('./data_files/GPL6244-tbl-1.txt', header=None, delimiter="\t")

# df_all_values.columns = ['id', 'gb_list', 'spot_id', 'seq_name', 'range_gb', 'range_strand', 'range_start',
#                          'range_stop',
#                          'total_probes', 'gene_assignment', 'mrna_assignment', 'category']

# df_all_values['gene_assignment'] = gene_assignment(df_all_values["gene_assignment"])

# data = join_column_dispersion(df_all_values, "GSM958226-tbl-1.txt")

# Removed the rows where disp == -1.0 , Null
# df_remove = data.loc[data.disp != -1.0]
