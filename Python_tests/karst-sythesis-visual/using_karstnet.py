import karstnet as kn
import pandas as pd
import numpy as np
from sklearn import preprocessing
from sklearn import cluster
from sklearn import decomposition
import networkx as nx
import matplotlib.pyplot as plt
from typing import Union, Any
import matplotlib as mpl


annotation : mpl.pyplot.Annotation = None


def annotatePlot(ax : plt.Axes,  text : str, coords : list):
    global annotation
    if annotation is not None :
        if annotation.axes != ax:
            annotation.remove()
            annotation = ax.annotate("", xy=(0,0), xytext=(-50,20), textcoords="offset pixels",
                            bbox=dict(boxstyle="round", fc="w"),
                            arrowprops=dict(arrowstyle="->"), zorder=100)
    else :
        annotation = ax.annotate("", xy=(0,0), xytext=(-50,20),textcoords="offset pixels",
                bbox=dict(boxstyle="round", fc="w"),
                arrowprops=dict(arrowstyle="->"), zorder=100)
    annotation.xy = coords
    annotation.set_text(text)
    # annotation.get_bbox_patch().set_facecolor(plt.cm.RdYlGn(plt.Normalize(1,4)(c[ind["ind"][0]])))
    annotation.get_bbox_patch().set_alpha(0.4)

def displayDataOnHover(event : mpl.backend_bases.Event, karsts_names : list, fig : plt.Figure):
    global annotation
    newAnnotationVis = False
    ax : mpl.pyplot.Axes = event.inaxes
    if ax is not None:
        contained, ind = ax.collections[0].contains(event)
        closestIndex = -1
        closestDist = 10000000
        best_coords = [-1, -1]
        for i in ind['ind']:
            x, y = ax.collections[0].get_offsets()[i]
            if closestDist > x**2 + y**2:
                closestDist = x**2 + y**2
                closestIndex = i
                best_coords = [x, y]
        if contained:
            annotatePlot(ax, karsts_names[closestIndex], best_coords)
            newAnnotationVis = True

        if annotation is not None:
            annotation.set_visible(newAnnotationVis)
        fig.canvas.draw_idle()

def getAllKarstsData(_desired_columns: dict) -> pd.DataFrame:
    filename = "data/karst_data_P_Collon.csv"
    df = pd.read_csv(filename, sep=";", decimal=',')
    # df.drop(columns=['s', 'n', 'L'], inplace=True)
    if isinstance(_desired_columns, dict):
        desired_columns = [col for col in _desired_columns if col in df.columns and _desired_columns[col]]
        df = df[desired_columns]
    return df


def normalizeData(df: pd.DataFrame, scaler_type: Any = preprocessing.MinMaxScaler) -> [pd.DataFrame, preprocessing.MinMaxScaler]:
    columns = df.columns
    x = df.drop(columns=['KarstName']).astype(float).values  # returns a numpy array
    min_max_scaler = scaler_type()
    x_scaled = min_max_scaler.fit_transform(x)
    df = pd.concat([df['KarstName'], pd.DataFrame(x_scaled)], axis=1)
    df.columns = columns
    return df, min_max_scaler


def decomposeData(df: pd.DataFrame, number_of_components: int = 2) -> [pd.DataFrame, decomposition.PCA]:
    columns = df.columns
    x = df.drop(columns=['KarstName']).astype(float).values  # returns a numpy array
    pca = decomposition.PCA(n_components=number_of_components)
    x_scaled = pca.fit_transform(x)
    df = pd.concat([df['KarstName'], pd.DataFrame(x_scaled)], axis=1)
    df.columns = ['KarstName'] + [f"Axis{i + 1}" for i in range(pca.n_components)]
    return df, pca


def karstDataToDf(karst: kn.KGraph, scaler: Union[preprocessing.MinMaxScaler, decomposition.PCA],
                  _desired_columns: dict) -> np.array:
    data_dict = karst.characterize_graph()

    # On the simplified graph
    nb_nodes_simplified = nx.number_of_nodes(karst.graph_simpl)
    nb_edges_simplified = nx.number_of_edges(karst.graph_simpl)
    nb_connected_components = nx.number_connected_components(
        karst.graph_simpl)

    nb_cycles = nb_edges_simplified - nb_nodes_simplified + nb_connected_components

    alpha = nb_cycles / (2 * (nb_nodes_simplified) - 5)
    beta = nb_edges_simplified / (nb_nodes_simplified)
    gamma = nb_edges_simplified / (3 * (nb_nodes_simplified - 2))
    connect_degree = (alpha/0.25 + (beta-1)/0.5 + (gamma - 0.33)/0.17)/3

    total_length = np.sum(karst.br_lengths) / 1000
    number_nodes = nx.number_of_nodes(karst.graph)
    number_edges = nx.number_of_edges(karst.graph)
    df = pd.DataFrame(np.array([total_length,
                                data_dict['orientation entropy'],
                                number_nodes,
                                number_edges,
                                data_dict['mean length'],
                                data_dict['length entropy'],
                                data_dict['cv length'],
                                data_dict['tortuosity'],
                                data_dict['mean degree'],
                                data_dict['cv degree'],
                                data_dict['correlation vertex degree'],
                                data_dict['aspl'],
                                data_dict['cpd'],
                                nb_nodes_simplified,
                                nb_edges_simplified,
                                nb_connected_components,
                                nb_cycles,
                                alpha,
                                beta,
                                gamma,
                                connect_degree,
                                None # The theta parameter is not computed yet...
                                ]).reshape(1, -1),
                      columns=['L', 'Ho', 'n', 's', 'len', 'Hlen', 'CVlen', 't', 'k', 'CVk', 'rk', 'SPL', 'CPD', 'N', 'S', 'p', 'Ncycl', 'alpha', 'beta', 'gamma', 'Dc', 'theta'])

    if isinstance(_desired_columns, dict):
        desired_columns = [col for col in _desired_columns if col in df.columns and _desired_columns[col]]
        df = df[desired_columns]
    columns = df.columns.insert(0, 'KarstName')
    x = df.astype(float).values  # returns a numpy array
    x_scaled = scaler.transform(x)
    df = pd.concat([pd.DataFrame([['myKarst']], columns=['KarstName']), pd.DataFrame(x_scaled)], axis=1)
    if isinstance(scaler, preprocessing.MinMaxScaler):
        df.columns = columns
    # data = scaler.transform(data.reshape(1, -1))
    # data = np.concatenate([["myKarst"], data[0]])
    data = df.to_numpy()[0]
    return data


def compareToAllKarsts(karst_data: np.array, all_karsts: pd.DataFrame) -> np.array:
    distances = {}
    for karst in all_karsts.iloc:
        distances[karst['KarstName']] = {
            'dist': np.linalg.norm(np.array(karst[1:].astype(float) - karst_data[1:].astype(float)))}
    distances = dict(sorted(distances.items(), key=lambda item: item[1]['dist']))
    return distances


def plotKmeans(all_karsts: pd.DataFrame, nb_centers: int, scaler: Union[preprocessing.MinMaxScaler, decomposition.PCA],
               karst_names: list, **kw) -> np.array:
    kmeans = cluster.KMeans(n_clusters=nb_centers)
    labels = kmeans.fit_predict(all_karsts)
    columns = all_karsts.columns
    if isinstance(scaler, preprocessing.MinMaxScaler):
        all_karsts = pd.DataFrame(scaler.inverse_transform(all_karsts), columns=columns)
    colors = labels
    linesizes = [0] * (len(colors) - 1) + [1]
    linecolors = (0, 0, 0, 0)
    alpha = [0.5] * (len(colors) - 1) + [1]
    fig, axes = plt.subplots(len(columns), len(columns), gridspec_kw={'wspace': 0, 'hspace': 0}, **kw)

    for i in range(len(columns)):
        for j in range(len(columns)):
            ax = axes[i][j]
            ax.scatter(all_karsts[columns[j]], all_karsts[columns[i]], c=colors, cmap=plt.cm.jet, linewidths=linesizes,
                       alpha=alpha, marker="o", s=15, edgecolors=linecolors)
            ax.get_xaxis().set_ticks([])
            ax.get_yaxis().set_ticks([])
            ax.set_xmargin(0)
            ax.set_ymargin(0)
            x_line = all_karsts[columns[j]].tolist()[-1]
            y_line = all_karsts[columns[i]].tolist()[-1]
            ax.axvline(x=x_line, c='red', linewidth=0.5, zorder=0, alpha=0.5)
            ax.axhline(y=y_line, c='red', linewidth=0.5, zorder=0, alpha=0.5)

            if i == len(columns) - 1:
                ax.set_xlabel(columns[j])
            if j == 0:
                ax.set_ylabel(columns[i])
    fig.canvas.mpl_connect("motion_notify_event", lambda evt: displayDataOnHover(evt, karst_names, fig))
    plt.tight_layout()
    plt.show()
    return labels


def plotDBScan(all_karsts: pd.DataFrame, scaler: Union[preprocessing.MinMaxScaler, decomposition.PCA],
               karst_names: list, **kw) -> np.array:
    dbscan = cluster.DBSCAN(eps=1.0, min_samples=3)
    labels = dbscan.fit_predict(all_karsts)
    columns = all_karsts.columns
    if isinstance(scaler, preprocessing.MinMaxScaler):
        all_karsts = pd.DataFrame(scaler.inverse_transform(all_karsts), columns=columns)
    colors = labels
    linesizes = [0] * (len(colors) - 1) + [1]
    linecolors = (0, 0, 0, 0)
    alpha = [0.5] * (len(colors) - 1) + [1]

    fig, axes = plt.subplots(len(columns), len(columns), gridspec_kw={'wspace': 0, 'hspace': 0}, **kw)
    for i in range(len(columns)):
        for j in range(len(columns)):
            ax = axes[i][j]
            ax.scatter(all_karsts[columns[j]], all_karsts[columns[i]], c=colors, cmap=plt.cm.jet, linewidths=linesizes,
                       alpha=alpha, marker="o", s=15, edgecolors=linecolors)
            ax.get_xaxis().set_ticks([])
            ax.get_yaxis().set_ticks([])
            ax.set_xmargin(0)
            ax.set_ymargin(0)
            ax.axhline(y=all_karsts[columns[i]].iloc[-1], c='red', linewidth=0.5, zorder=0, alpha=0.5)
            ax.axvline(x=all_karsts[columns[j]].iloc[-1], c='red', linewidth=0.5, zorder=0, alpha=0.5)

            if i == len(columns) - 1:
                ax.set_xlabel(columns[j])
            if j == 0:
                ax.set_ylabel(columns[i])
    fig.canvas.mpl_connect("motion_notify_event", lambda evt: displayDataOnHover(evt, karst_names, fig))
    plt.tight_layout()
    plt.show()
    return labels


def plotDBScanUntilOutlierExists(all_karsts: pd.DataFrame, scaler: Union[preprocessing.MinMaxScaler, decomposition.PCA],
                                 karst_names: list, **kw) -> np.array:
    eps = 50.0
    nb_labels = 1
    labels = []
    while nb_labels == 1 and eps > 0.001:
        dbscan = cluster.DBSCAN(eps=eps, min_samples=3)
        labels = dbscan.fit_predict(all_karsts)
        nb_labels = len(np.unique(labels))
        if nb_labels == 1:
            eps *= 0.9
    columns = all_karsts.columns
    if isinstance(scaler, preprocessing.MinMaxScaler):
        all_karsts = pd.DataFrame(scaler.inverse_transform(all_karsts), columns=columns)
    colors = labels
    linesizes = [0] * (len(colors) - 1) + [1]
    linecolors = (0, 0, 0, 0)
    alpha = [0.5] * (len(colors) - 1) + [1]

    fig, axes = plt.subplots(len(columns), len(columns), gridspec_kw={'wspace': 0, 'hspace': 0}, **kw)
    for i in range(len(columns)):
        for j in range(len(columns)):
            ax = axes[i][j]
            ax.scatter(all_karsts[columns[j]], all_karsts[columns[i]], c=colors, cmap=plt.cm.jet, linewidths=linesizes,
                       alpha=alpha, marker="o", s=15, edgecolors=linecolors)
            ax.get_xaxis().set_ticks([])
            ax.get_yaxis().set_ticks([])
            ax.set_xmargin(0)
            ax.set_ymargin(0)
            ax.axhline(y=all_karsts[columns[i]].iloc[-1], c='red', linewidth=0.5, zorder=0, alpha=0.5)
            ax.axvline(x=all_karsts[columns[j]].iloc[-1], c='red', linewidth=0.5, zorder=0, alpha=0.5)

            if i == len(columns) - 1:
                ax.set_xlabel(columns[j])
            if j == 0:
                ax.set_ylabel(columns[i])
    fig.canvas.mpl_connect("motion_notify_event", lambda evt: displayDataOnHover(evt, karst_names, fig))
    print("DBScan with eps = " + str(eps) + " -- " + str(nb_labels) + " classes")
    plt.tight_layout()
    plt.show()
    return labels


def plotAgglomerativeCluster(all_karsts: pd.DataFrame, nb_centers: int,
                             scaler: Union[preprocessing.MinMaxScaler, decomposition.PCA], karst_names: list, **kw) -> np.array:
    model = cluster.AgglomerativeClustering(n_clusters=nb_centers, linkage="complete")
    labels = model.fit_predict(all_karsts)
    columns = all_karsts.columns
    if isinstance(scaler, preprocessing.MinMaxScaler):
        all_karsts = pd.DataFrame(scaler.inverse_transform(all_karsts), columns=columns)
    colors = labels
    linesizes = [0] * (len(colors) - 1) + [1]
    linecolors = (0, 0, 0, 0)
    alpha = [0.5] * (len(colors) - 1) + [1]

    fig, axes = plt.subplots(len(columns), len(columns), gridspec_kw={'wspace': 0, 'hspace': 0}, **kw)
    for i in range(len(columns)):
        for j in range(len(columns)):
            ax = axes[i][j]
            ax.scatter(all_karsts[columns[j]], all_karsts[columns[i]], c=colors, cmap=plt.cm.jet, linewidths=linesizes,
                       alpha=alpha, marker="o", s=15, edgecolors=linecolors)
            ax.get_xaxis().set_ticks([])
            ax.get_yaxis().set_ticks([])
            ax.set_xmargin(0)
            ax.set_ymargin(0)
            ax.axhline(y=all_karsts[columns[i]].iloc[-1], c='red', linewidth=0.5, zorder=0, alpha=0.5)
            ax.axvline(x=all_karsts[columns[j]].iloc[-1], c='red', linewidth=0.5, zorder=0, alpha=0.5)

            if i == len(columns) - 1:
                ax.set_xlabel(columns[j])
            if j == 0:
                ax.set_ylabel(columns[i])
    fig.canvas.mpl_connect("motion_notify_event", lambda evt: displayDataOnHover(evt, karst_names, fig))
    plt.tight_layout()
    plt.show()
    return labels


def main():
    # filename = "data/karst"
    filename = "data/Huttes"
    # filename = "data/Sakany"
    # filename = "data/superimposed"
    # filename = "data/spongework"
    # filename = "data/rectilinear"
    # filename = "data/gorge"
    karst = kn.from_nodlink_dat(filename)
    # karst.plot3()

    desired_columns = {'KarstName': 1,
                       's': 0,
                       'n': 0,
                       'Ho': 1,
                       'len': 1,
                       'Hlen': 1,
                       'CVlen': 1,
                       't': 1,
                       'k': 1,
                       'CVk': 1,
                       'rk': 1,
                       'SPL': 1,
                       'CPD': 1,
                       'N' : 1,
                       'S' : 1,
                       'p' : 1,
                       'Ncycl' : 1,
                       'alpha' : 1,
                       'beta' : 1,
                       'gamma' : 1,
                       'Dc' : 1,
                       'theta' : 0 # Cannot be set for now
                       }

    df, scaler = normalizeData(getAllKarstsData(desired_columns), preprocessing.QuantileTransformer)
    df_pca, pca = decomposeData(getAllKarstsData(desired_columns))
    karst_data = karstDataToDf(karst, scaler, desired_columns)
    karst_data_pca = karstDataToDf(karst, pca, desired_columns)
    print(compareToAllKarsts(karst_data, df))
    print(compareToAllKarsts(karst_data_pca, df_pca))

    all_karsts = pd.concat([df, pd.DataFrame(karst_data.reshape(1, -1), columns=list(df))], axis=0, ignore_index=True)
    names = all_karsts["KarstName"]
    all_karsts.drop(columns=["KarstName"], inplace=True)

    all_karsts_pca = pd.concat([df_pca, pd.DataFrame(karst_data_pca.reshape(1, -1), columns=list(df_pca))], axis=0,
                               ignore_index=True) \
        .drop(columns=["KarstName"]).astype(float)

    kmeans_labels = plotKmeans(all_karsts, 3, scaler, names.to_list())
    myIndex = kmeans_labels[-1]
    print("Using Kmeans, our karst looks like : ", names[kmeans_labels == myIndex].tolist())
    """
        dbscan_labels = plotDBScanUntilOutlierExists(all_karsts, scaler, names.to_list())
        myIndex = dbscan_labels[-1]
        print("Using DBScan, our karst looks like : ",names[dbscan_labels == myIndex].tolist())
    
        agglomerative_labels = plotAgglomerativeCluster(all_karsts, 3, scaler, names.to_list())
        myIndex = agglomerative_labels[-1]
        print("Using AgglomerativeClustering, our karst looks like : ",names[agglomerative_labels == myIndex].tolist())"""

    # PCA part
    kmeans_labels = plotKmeans(all_karsts_pca, 3, pca, names.to_list())
    myIndex = kmeans_labels[-1]
    print("Using Kmeans, our karst looks like : ", names[kmeans_labels == myIndex].tolist())
    """
    dbscan_labels = plotDBScanUntilOutlierExists(all_karsts_pca, pca, names.to_list())
    myIndex = dbscan_labels[-1]
    print("Using DBScan, our karst looks like : ",names[dbscan_labels == myIndex].tolist())

    agglomerative_labels = plotAgglomerativeCluster(all_karsts_pca, 3, pca, names.to_list())
    myIndex = agglomerative_labels[-1]
    print("Using AgglomerativeClustering, our karst looks like : ", names[agglomerative_labels == myIndex].tolist())"""


if __name__ == '__main__':
    main()
