

import numpy as np
import pandas as pd

import PySUSSIX.f90sussix.f90sussix as f90sussix
import PySUSSIX.crossref as crossref



def datspe(x,px,y,py,zeta,pzeta,number_of_harmonics = 5,method = 'hanning',return_values = False):
    kwargs = {  'fortran_unit_number'       : 90,
                'dimension_of_phase_space'  : 3,
                'flag_for_real_signal'      : 1,
                'initial_turn_number'       : 1,
                'final_turn_number'         : len(x),
                'number_of_turns'           : len(x),
                'flag_for_windowing'        : {'hanning':1,'rectangular':2}[method],
                'number_of_harmonics'       : number_of_harmonics,
                'flag_for_full_analysis'    : 1}
    

    fortran_kwargs = {key:kwargs[crossref.fort2py[key]] for key in ['iunit','idam','ir','nt1','nt2','nturn','imeth','narm','iana']}
    fortran_kwargs['x']     = x
    fortran_kwargs['xp']    = px
    fortran_kwargs['y']     = y
    fortran_kwargs['yp']    = py
    fortran_kwargs['s']     = zeta
    fortran_kwargs['sp']    = pzeta
    fortran_kwargs['n_points'] = len(x)

    # run subroutine
    f90sussix.datspe(*fortran_kwargs.values())

    # The values are stored in the COMMONS (fortran database)
    if return_values:
        return pd.DataFrame({   'zxpes':f90sussix.fcoe.zxpes[:number_of_harmonics],
                                'txa'  :f90sussix.tune.txa[:number_of_harmonics],
                                'zypes':f90sussix.fcoe.zxpes[:number_of_harmonics],
                                'tya'  :f90sussix.tune.tya[:number_of_harmonics],
                                'zspes':f90sussix.fcoe.zxpes[:number_of_harmonics],
                                'tsa'  :f90sussix.tune.tya[:number_of_harmonics]})
                        

