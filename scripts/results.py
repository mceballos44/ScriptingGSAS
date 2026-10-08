# Functions here will extract parameters from sequential refinement
# Also set up nice tables
# We want the data to be setup like this
# Sample, T, A, sigmaA, C, sigmaC, Diso, Dsigma, strain,sigmaS, Rwp
import os
import G2script as G2sc
from pathlib import Path
import pandas as pd
from datetime import datetime
# Import directories
from scripts.config import setup,BASE_DIR,DATA_DIR,PROJECT_DIR,OUTPUT_DIR,CONTROLS_DIR
from scripts.config import BACKGROUND_FILE,CONTROLS_FILE,MASK_FILE,INSTRUMENT_FILE
# Will always use this cif file for the phase
ice_cif = DATA_DIR / "ice.cif"

columns=['Sample','T','A','sigmaA','C','sigmaC',
                 'D11','sigma_D11','D33','sigma_D33',
                 'Astrain','sigma_Astrain','Cstrain',
                 'sigma_Cstrain','Rwp']

def extract_data(gpx,sample_name):
    """
    This function should extract the data from the current gpx and 
    store data in dataframe
    """
    rows = []
    seq = gpx.seqref()

    # print(f'\n\nParamlist: {seq.get_ParmList(0)}')
    # print(f'\n\nCell and ESD format: {seq.get_cell_and_esd('ice',0)}')
    # We can use dictionaries instead of constantly appending data
    for x,_ in enumerate(seq.histograms()):
        cell, cellESD, _ = seq.get_cell_and_esd('ice',x)
        ref_data = seq.RefData(x)
        seq_results = ref_data[0]
        row = {
            'Sample': sample_name,
            'T': seq_results['parmDict'][f':{x}:Temperature'],
            'A': cell[0],
            'sigmaA': cellESD[0],
            'C': cell[2],
            'sigmaC': cellESD[2],
            'D11': seq.get_Variable(x,f'0:{x}:D11')[0],
            'sigma_D11': seq.get_Variable(x,f'0:{x}:D11')[1],
            'D33': seq.get_Variable(x,f'0:{x}:D33')[0],
            'sigma_D33': seq.get_Variable(x,f'0:{x}:D33')[1],
            'Mustrain_a': seq.get_Variable(x,f'0:{x}:Mustrain;a')[0],
            'sigma_Mustrain_a': seq.get_Variable(x,f'0:{x}:Mustrain;a')[1],
            'Mustrain_i': seq.get_Variable(x,f'0:{x}:Mustrain;i')[0],
            'sigma_Mustrain_i': seq.get_Variable(x,f'0:{x}:Mustrain;i')[1],
            'Rwp': ref_data[0]['Rvals']['Rwp']
        }
        rows.append(row)
    df = pd.DataFrame(rows)
    return df



def combine_data(sample_list):
    full_df = pd.concat(
        sample_list,
        ignore_index=True
    )
    
    # Save dataframe as readable csv/excel file
    return full_df

def save_data(df):
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M")
    excel_file = OUTPUT_DIR / f'seq_results_{timestamp}.xlsx'
    csv_file = OUTPUT_DIR / 'seq_results.csv'
    
    df.to_excel(
        excel_file,
        index=False
    )
    df.to_csv(
        csv_file,
        index=False
    )
    print(f"Saved results to:")
    print(excel_file)
    print(csv_file)
    return excel_file
# import scripts.config as sp
# import image_processing as ip
# import refinement as rf
# gpx = sp.setup()
# ip.integrate_images(gpx=gpx,samples ='AFP/AFP3')
# ip.remove_unfrozen(gpx=gpx)
# ip.remove_orphan_images(gpx=gpx)

# ip.assign_phase_one(gpx=gpx,phase_name='ice')

# rf.initial_refine(gpx=gpx)
# rf.first_seq_refine(gpx=gpx)
# rf.additional_seq_refine(gpx=gpx)
# extract_data(gpx=gpx,sample_name='AFP3')
