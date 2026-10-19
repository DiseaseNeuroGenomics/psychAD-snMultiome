import sys
import os

import pandas as pd
import numpy as np
import mofax as mfx
from mofapy2.run.entry_point import entry_point
from pathlib import Path

os.chdir('/sc/arion/projects/CommonMind/roussp01a/snmulti/DE/files/donor/new')
output_dir='MOFA'
Path(output_dir).mkdir(parents=True, exist_ok=True)

Sel=pd.read_csv('sel.csv')

chr_block=pd.read_csv('chr_block.csv')
Blubs={cell:pd.read_csv(f'{cell}_RandEf.csv.gz',index_col=0) for cell in Sel.cell.unique()}
for nf in [10]:
    for threshold in ['top25']:
        sel=Sel.loc[Sel.threshold==threshold]
        for rd in chr_block['round'].unique():
            rms=[f'chr{x}' for x in chr_block.loc[chr_block['round'].eq(rd)*chr_block.threshold.eq(threshold)].chr.values]
            sel1=sel.loc[~sel.chr.isin(rms)]
            df=list()
            for cell in sel1.cell.unique():
                gs=sel1.loc[sel1.cell.eq(cell)].gene.values
                z=Blubs[cell].loc[gs].reset_index(names='feature').melt(id_vars='feature', var_name='sample')
                z['view']=cell
                z['feature']=cell+':'+z['feature']
                df.append(z)
            data_dt=pd.concat(df)
            ent = entry_point()
            ent.set_data_df(data_dt, 
                    likelihoods = ["gaussian"]*data_dt.view.unique().shape[0])
            ent.set_data_options(scale_views = True)
            ent.set_model_options(factors = nf, spikeslab_weights = True, 
                          ard_weights = True)
            ent.set_train_options(convergence_mode = "slow", 
                          dropR2 = -1, gpu_mode = False, seed = 1)

            ent.build()
            ent.run()
            ent.save(outfile=f'{output_dir}/MOFA_{threshold}_{rd}_{nf}.hdf5') 

