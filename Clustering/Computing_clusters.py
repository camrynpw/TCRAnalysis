#see also: def compute_leiden_clusters_from_vecs(vecs, num_nbrs=5, random_seed=20):
# in phil_hacking.py
#
import os
import igraph as ig 
import numpy as np
import pandas as pd
import faiss
from timeit import default_timer as timer
import leidenalg
import random
import umap

# read in SNE lists, make list of unique TCRs
os.chdir('filepath/immune_response_detection')
from raptcr.hashing import TCRDistEncoder

# Read in the sequences file you want to cluster
df = pd.read_csv('filepath/sequences.csv')
# Split the "Sequence" column into "v_call" and "junction_aa" columns
df[['v_call', 'junction_aa']] = df['Sequence'].str.split('_', n=1, expand=True)

# encode the sequences from the df
encoder = TCRDistEncoder()
encoder.fit()
vecs = encoder.transform(df) # encoded vectors

# Save the embeddings into two dimensions with umap for later plotting
# adjust neighbors and min_dist based on data
umap_embeddings = umap.UMAP(n_neighbors=35, min_dist=0.6, metric='euclidean').fit_transform(vecs)
df['umap_1'] = umap_embeddings[:,0]
df['umap_2'] = umap_embeddings[:,1]

# create the graph
g = ig.Graph(directed=False) # undirected graph
num_tcrs = vecs.shape[0]
g.add_vertices(num_tcrs)

# now we use faiss to quickly find pairs of tcrs that are within a certain distance
idx = faiss.IndexFlatL2(vecs.shape[1])
idx.add(vecs)

radius = 12.5 # test to see what is good clustering radius
min_group_size = 2 #adjust this - what should it be?
start = timer()

# lims is an array that represents the neighbors/close ones with their indcs?
# basically a way of representing which ones are close together
lims, D, I = idx.range_search(vecs, radius)
print(f'range_search took {timer()-start:.2f}')
print('total nbrs:', radius, I.shape, flush=True)

# should be the length of the number of tcrs + 1 logically
assert lims.shape[0] == num_tcrs + 1
# calculating the neighbors for use in the graph
nnbrs = lims[1:]-lims[:-1]
I0 = np.repeat(np.arange(num_tcrs), nnbrs.astype(int))
start = timer()
assert I0.shape == I.shape
g.add_edges(list(zip(I0,I)))

print(f'add nbrs took {timer()-start:.2f}', flush=True)

# Leiden
random_seed = 20
partition_type = leidenalg.RBConfigurationVertexPartition
part = leidenalg.find_partition(g, partition_type, seed=random_seed)
num_vertices = g.vcount()
leiden_clusters = np.array(part.membership)

# Initialize the 'leiden' column in the DataFrame
df['leiden'] = leiden_clusters


# iteratively find the highest-degee vertex and remove it and its nbrs
# these are the groups
NOGROUP = -1
group_numbers = np.full((num_tcrs,), NOGROUP, dtype=int)
vertex2tcr = list(range(num_tcrs)) # map from graph vertices to tcr numbers
centers = []
while True:
    assert len(vertex2tcr) == g.vcount()
    degreelist = g.degree()
    center = np.argmax(degreelist)
    nbrs = sorted(set(g.neighbors(center)))
    assert center in nbrs
    print(g.vcount(), center, degreelist[center], g.maxdegree(), len(nbrs))
    tcr_nbrs = [vertex2tcr[x] for x in nbrs]
    tcr_center = vertex2tcr[center]
    assert np.all(group_numbers[tcr_nbrs]==NOGROUP)
    group_numbers[tcr_nbrs] = len(centers)
    centers.append(tcr_center)
    g.delete_vertices(nbrs)
    for v in sorted(nbrs, reverse=True):
        vertex2tcr.pop(v)
    if g.vcount()==0:
        break
    if g.maxdegree()+1 < min_group_size:
        break

mask = group_numbers==NOGROUP
if mask.sum():
    next_group = max(group_numbers)+1
    group_numbers[mask] = np.arange(next_group, next_group+mask.sum())

df['simple_cluster'] = group_numbers

# centers is the list of cluster centers
# group_numbers is a big numpy array of the cluster group membership

#Save the results in a file
outfile = 'filepath/Sequence_clusters.csv'
df.to_csv(outfile, sep=',', index=False)
print('made:', outfile)