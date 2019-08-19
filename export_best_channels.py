import numpy as np
import pandas as pd
from phycontrib.template import TemplateController
from phy.utils._misc import _read_python
from phy.gui.qt import create_app
import sys
from os import chdir

## REMEMBER PHY **HATES** EXTRA CSV FILES


# path = r'G:\My Drive\Documents\UMass\Healey Lab\AHNMDA01\L Hemi 180606\Playback recordings\alldata\\'
# output_prename = "AHSKF01_R_kilosort"
def export_best_channels(path):
    create_app()

    chdir(path)

    all_clusters = np.load(path + '\\spike_clusters.npy')
    cluster_quality = pd.read_csv(path + '\\cluster_groups.csv', sep = "\t")

    tc = TemplateController(**_read_python(path + '\\params.py'))

    cluster_best_dict = dict()
    for cluster in all_clusters:
        cluster_best_dict.update({cluster: tc.get_best_channel(cluster)})

    cluster_best_df = pd.DataFrame.from_dict(cluster_best_dict, orient='index')
    cluster_best_df.reset_index(level=0, inplace=True)
    cluster_best_df.columns = ['Cluster_id', 'Best_channel']
    cluster_best_df.to_csv(path + '\\best_channels.csv', index=False)


if __name__ == '__main__':
    path = str(sys.argv[1])
    export_best_channels(path)