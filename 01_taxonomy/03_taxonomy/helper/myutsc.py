import numpy as np
import pandas as pd
from anndata import AnnData
from typing import Optional
import h5py
from scipy import sparse
import scanpy as sc

def injest_from_same_cohort(ad_ref,ad_query,use_rep='X_pca_harmony',n_pcs=10,
                            epsilon=0.1,obs=None):
    # adapted from scanpy injest scanpy.tl.ingest
    from scipy.sparse import issparse
    from sklearn.utils import check_random_state
    from pynndescent import NNDescent
    from pynndescent.rp_trees import make_forest
    from umap.umap_ import INT32_MAX, INT32_MIN

    neighbors_key = "neighbors"
    neighbors = sc._utils.NeighborsView(ad_ref, neighbors_key)
    n_neighbors = neighbors["params"]["n_neighbors"]
    rep=ad_ref.obsm[use_rep][:,:n_pcs]
    if "metric_kwds" in neighbors["params"]:
        metric_kwds = neighbors["params"]["metric_kwds"]
    else:
        metric_kwds = {}

    metric = neighbors["params"]["metric"]
    neigh_random_state = neighbors["params"].get("random_state", 0)
    distances=neighbors["distances"]
    first_col = np.arange(distances.shape[0])[:, None]
    init_indices = np.hstack((first_col, np.stack(distances.tolil().rows)))
    nnd_idx = NNDescent(
            data=rep,
            metric=metric,
            metric_kwds=metric_kwds,
            n_neighbors=n_neighbors,
            init_graph=init_indices,
            random_state=neigh_random_state)
    current_random_state = check_random_state(nnd_idx.random_state)
    nnd_idx._rp_forest = make_forest(
            nnd_idx._raw_data,
            nnd_idx.n_neighbors,
            nnd_idx.n_search_trees,
            nnd_idx.leaf_size,
            nnd_idx.rng_state,
            current_random_state,
            nnd_idx.n_jobs,
            nnd_idx._angular_trees,
        )
    random_state=0
    random_state = check_random_state(random_state)
    rng_state = random_state.randint(INT32_MIN, INT32_MAX, 3).astype(np.int64)

    test = ad_query.obsm[use_rep][:,:n_pcs]

    nnd_idx.search_rng_state = rng_state
    _indices, _distances = nnd_idx.query(test, n_neighbors, epsilon)
    ad_query.obs[f'pred_{obs}']=[ad_ref.obs[obs].iloc[inds].mode()[0] for inds in _indices]
    return ad_query
brain_cell_markers={"GLU":["RBFOX3", "MEG3", "SLC17A7"],
"GABA":["RBFOX3", "MEG3", "GAD1", "GAD2", "GRIK1"],
"Astro":["SLC1A3", "GFAP", "APOE", "SLC1A2", "SLC14A1"],
"Oligo":["PLP1", "MAG", "MBP"],
"OPC": ["PDGFRA", "VCAN"],
"Endo":["FLT1", "CLDN5"],
"Mural":["PDGFRB"],
"Microg": ["TGFBR1", "DOCK8", "CD74", "CSF1R", "MS4A6A", "PLXDC2"],
"Tcell":["CD3D", "CD3E", "CD3G","TRAC"],
"Erythroid":["GYPA","TFRC","ITGA4","HBB","HBA1"]}

max_268 = ["#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
        "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
        "#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
        "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
        "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
        "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
        "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
        "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
        "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
        "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
        "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
        "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
        "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C",
        "#83AB58", "#001C1E", "#D1F7CE", "#004B28", "#C8D0F6", "#A3A489", "#806C66", "#222800",
        "#BF5650", "#E83000", "#66796D", "#DA007C", "#FF1A59", "#8ADBB4", "#1E0200", "#5B4E51",
        "#C895C5", "#320033", "#FF6832", "#66E1D3", "#CFCDAC", "#D0AC94", "#7ED379", "#012C58",
        "#7A7BFF", "#D68E01", "#353339", "#78AFA1", "#FEB2C6", "#75797C", "#837393", "#943A4D",
        "#B5F4FF", "#D2DCD5", "#9556BD", "#6A714A", "#001325", "#02525F", "#0AA3F7", "#E98176",
        "#DBD5DD", "#5EBCD1", "#3D4F44", "#7E6405", "#02684E", "#962B75", "#8D8546", "#9695C5",
        "#E773CE", "#D86A78", "#3E89BE", "#CA834E", "#518A87", "#5B113C", "#55813B", "#E704C4",
        "#00005F", "#A97399", "#4B8160", "#59738A", "#FF5DA7", "#F7C9BF", "#643127", "#513A01",
        "#6B94AA", "#51A058", "#A45B02", "#1D1702", "#E20027", "#E7AB63", "#4C6001", "#9C6966",
        "#64547B", "#97979E", "#006A66", "#391406", "#F4D749", "#0045D2", "#006C31", "#DDB6D0",
        "#7C6571", "#9FB2A4", "#00D891", "#15A08A", "#BC65E9", "#FFFFFE", "#C6DC99", "#203B3C",
        "#671190", "#6B3A64", "#F5E1FF", "#FFA0F2", "#CCAA35", "#374527", "#8BB400", "#797868",
        "#C6005A", "#3B000A", "#C86240", "#29607C", "#402334", "#7D5A44", "#CCB87C", "#B88183",
        "#AA5199", "#B5D6C3", "#A38469", "#9F94F0", "#A74571", "#B894A6", "#71BB8C", "#00B433",
        "#789EC9", "#6D80BA", "#953F00", "#5EFF03", "#E4FFFC", "#1BE177", "#BCB1E5", "#76912F",
        "#003109", "#0060CD", "#D20096", "#895563", "#29201D", "#5B3213", "#A76F42", "#89412E",
        "#1A3A2A", "#494B5A", "#A88C85", "#F4ABAA", "#A3F3AB", "#00C6C8", "#EA8B66", "#958A9F",
        "#BDC9D2", "#9FA064", "#BE4700", "#658188", "#83A485", "#453C23", "#47675D", "#3A3F00",
        "#061203", "#DFFB71", "#868E7E", "#98D058", "#6C8F7D", "#D7BFC2", "#3C3E6E", "#D83D66",
        "#2F5D9B", "#6C5E46", "#D25B88", "#5B656C", "#00B57F", "#545C46", "#866097", "#365D25",
        "#252F99", "#00CCFF", "#674E60", "#FC009C", "#92896B"]


