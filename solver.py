import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
import sys
import os

if __name__ == "__main__":

    nproc = len(os.listdir('OUT/'))
    nsamples = 1021

    nfeat = -1
    nc = -1
    nloc = -1
    myloc = -1

    fdata = None
    counter = None
    centers = None
    clustId = np.zeros(nsamples)

    file_string = "./OUT/FinalOutId"

    for p_id in range(nproc):
        current_file = file_string + str(p_id)
        f = open(current_file)
        f.readline()
        params = f.readline().split()
        myloc = int(params[3])
        if p_id == 0:

            nc = int(params[1][:-1])
            nfeat = int(params[5][:-1])
            nloc = myloc

            fdata = np.zeros((nsamples, nfeat))
            counter = np.zeros(nc)
            centers = np.zeros((nc, nfeat))


            f.readline()
            for line_idx in range(nc):
                counts = f.readline().split()
                counter[line_idx] = int(counts[3])

        for line_idx in range(myloc):
            clustId[p_id * nloc + line_idx] = int(f.readline())

        f.readline()

        for line_idx in range(myloc):
            data_pt = f.readline().split()
            data_pt = list(map(float, data_pt))
            fdata[p_id * nloc + line_idx] = np.array(data_pt)

        f.close()



    for i in range(nc):
        cluster_data = fdata[clustId == i]
        centers[i] = np.mean(cluster_data, axis=0)


    tol = 1e-10
    n_init=1
    kmeans = KMeans(n_clusters=nc,tol=tol, n_init=n_init, init=centers)
    kmeans.fit(fdata)

    correct_labels = np.sum(clustId == kmeans.labels_)
    print("number of correct labels = ",correct_labels, "/", nsamples)
