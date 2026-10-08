import os
import G2script as G2sc
from pathlib import Path
# Import directories
from scripts.config import setup,BASE_DIR,DATA_DIR,PROJECT_DIR,OUTPUT_DIR,CONTROLS_DIR
from scripts.config import BACKGROUND_FILE,CONTROLS_FILE,MASK_FILE,INSTRUMENT_FILE
import scripts.image_processing as ip
# Will always use this cif file for the phase
ice_cif = DATA_DIR / "ice.cif"

def initial_refine(gpx):
    """
    Refine a single pwdr set to a decent point before copying all
    data to other pwdr
    """

    hist = gpx.histogram(0)

    reflist = [
        {
            'set': {
                'Sample Parameters': ['DisplaceX'],
                'Background': {
                    'no. coeffs': 5,
                    'refine': True
                }
            }
        },

        {
            'set': {
                'Cell': True,
                'Mustrain': {
                    'type': 'uniaxial',
                    'refine': 'axial',
                    'direction': [0,0,1]
                    # C-axis has lower CTE generally, so we expect axial strain to be
                    # unique
                },
                'Mustrain':{
                    'type': 'uniaxial',
                    'refine': 'equatorial'
                }  
            }#,
            # 'once':{
            #     'Size':{
            #         'refine': True
            #     }
            # }
        },
        {
            'set': {
                'Cell': False,
                'Mustrain': {
                    'refine': False
                },
                'Background': {
                    'no. coeffs': 5,
                    'refine': False
                }
            },
            'clear': {
                'Sample Parameters': ['DisplaceX'],
            }
        },
        
        {
            'once': {
                'Instrument Parameters': ['Zero']
            }
        },
        
        {
            'set': {
                'Sample Parameters': ['DisplaceX'],
                'Background': {
                    'no. coeffs': 5,
                    'refine': True
                }
            },
            'set': {
                'Cell': True,
                'Mustrain': {
                    'refine': True
                }
            }
        },
        # Clear all parameters for sequential start
        {
            'set': {
                'Cell': False,
                'Mustrain': {
                    'refine': False
                },
                'Background': {
                    'refine': False
                }
            },
            'clear': {
                'Sample Parameters': ['DisplaceX']
            }
        }
    ]


    gpx.do_refinements(
        reflist,
        histogram=hist
    )
    gpx.save()
    return gpx

def copy_displaceX(gpx):
    """
    Copy DisplaceX from the first PWDR histogram
    to all remaining PWDR histograms.
    """

    # Get only PWDR histograms
    pwdr_list = [
        hist for hist in gpx.histograms()
        if hist.name.startswith("PWDR ")
    ]

    # First PWDR is our reference
    reference = pwdr_list[0]

    # Get DisplaceX value
    displace_x = reference.data["Sample Parameters"]["DisplaceX"]

    print(f"Reference: {reference.name}")
    print(f"DisplaceX: {displace_x}")

    # Copy to all remaining PWDR histograms
    for hist in pwdr_list[1:]:

        hist.data["Sample Parameters"]["DisplaceX"] = displace_x

        print(f"Copied DisplaceX -> {hist.name}")

    gpx.save()

    return gpx
def first_seq_refine(gpx):
    """
    Sets up sequential refinement and performs single sequential refinement step
    
    """
    
    ip.assign_phase_all(gpx,'ice') # Assign phase to all histograms
    # Need to copy all flags and parameters to all pwdrs and phase data
    gpx.copyHistParms(0,'all',['b','i','L'])
    phase = gpx.phases()[0]
    phase.copyHAPvalues(0,targethistlist='all')
    # Copy x displacement in sample params
    copy_displaceX(gpx=gpx)
    
    refs = {
            'set': {
                'HStrain': {
                    'refine': True
                }
            }
        }
    # Set Dij terms to refine
    gpx.set_refinement(refs)
    # Set sequential refinement controls and perform refinement
    gpx.set_Controls('sequential',gpx.histograms())
    gpx.set_Controls('cycles',10)
    gpx.set_Controls('seqCopy',True)
    
    # Set crystallite max and min size to reasonable value
    gpx.set_Controls(
        'parmMin',
        variable='0:*:Size;i',
        value=0.001
    )
    gpx.set_Controls(
    'parmMax',
    variable='0:*:Size;i',
    value=5.0
    )
    
    gpx.refine()
    gpx.save()
    return gpx

def additional_seq_refine(gpx):
    """
    Here we perform additional sequential refinement steps with more refined terms
    Goal here is to get better measures of lattice params and eventually size 
    """
    refs = [
        {
            'set':{
                'Sample Parameters': ['DisplaceX'],
                'Background': {
                    'refine': True
                }
            }
        },
        {
            'set':{
                'Mustrain':{
                    'refine': True
                }#,
                # 'Size':{
                #     'refine': True
                # }
            }
        }
    ]

    gpx.do_refinements(refs)
    
    #gpx.refine()
    gpx.save()
    return gpx
    
#Testing
# import scripts.config as sp
# import image_processing as ip
# gpx = sp.setup()
# ip.integrate_images(gpx=gpx,samples ='AFP/AFP3')
# ip.remove_unfrozen(gpx=gpx)
# ip.remove_orphan_images(gpx=gpx)

# ip.assign_phase_one(gpx=gpx,phase_name='ice')
# #ip.assign_phase_all(gpx=gpx,phase_name='ice')

# initial_refine(gpx=gpx)
# first_seq_refine(gpx=gpx)
# additional_seq_refine(gpx=gpx)
# print(gpx.get_Controls('sequential'))

# results = gpx.get_LastFitResults()

# for var in results:
#     if "Size" in str(var):
#         print(var)
#print(gpx.get_Frozen())