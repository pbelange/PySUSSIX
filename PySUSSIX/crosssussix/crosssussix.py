

import numpy as np
import pandas as pd

import PySUSSIX.f90sussix.f90sussix as f90sussix
import PySUSSIX.crosssussix.crossref as crossref



def datspe(x,px,y,py,zeta,pzeta,number_of_harmonics = 5,method = 'hanning',return_values = False):
    """
    THIS PROGRAM CALCULATES THE SPECTRUM OF A TRACKING DATA SET,
    THE TUNE AND THE LINES ARE CALCULATED WITH THE ROUTINE SPECTRUM
    """
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
                        


def ordres(nturns,tune_estimate = [0.31,0.32,0.0018],number_of_harmonics = 5,return_values = False):
    """
    Returning frequencies with jlkm indices where Q_eff = j*Q_x + k*Q_y + l*Q_z + m
    """
    assert not np.all(f90sussix.tune.txa == 0), 'No tune found, please run datspe first.'

    fortran_args = ['eps','narm','nrc','idam','iunit','nturn','tunex','tuney','tunez','istune','etune','amplitude','phase','ox','ax','oy','ay','os','as']

    kwargs = {  'tolerance_on_identification_of_frequencies'    :1e-10, 
                'number_of_harmonics'                           :number_of_harmonics,
                'maximum_order_of_combination_of_frequencies'   :10,
                'dimension_of_phase_space'                      :3,
                'fortran_unit_number'                           :90,
                'number_of_turns'                               :nturns,
                'tunex'                                         :tune_estimate[0],
                'tuney'                                         :tune_estimate[1],
                'tunez'                                         :tune_estimate[2],
                'flag_for_fundamental_frequencies'              :1,
                'allowed_distance_to_fundamental_frequencies'   :[1e-2,1e-2,1e-2],
                'amplitude_values'                              :np.zeros(14),
                'phase_values'                                  :np.zeros(14),
                'ox'                                            :np.zeros(600),
                'ax'                                            :np.zeros(600),
                'oy'                                            :np.zeros(600),
                'ay'                                            :np.zeros(600),
                'os'                                            :np.zeros(600),
                'as'                                            :np.zeros(600)}


    # Cross-referencing
    fortran_kwargs = {key:kwargs[crossref.fort2py[key]] for key in fortran_args}

    # run subroutine
    f90sussix.ordres(*fortran_kwargs.values())


    # The values are stored in ox,ax,amplitude,phase vectors
    # Returning frequencies with jlkm indices where Q_eff = j*Q_x + k*Q_y + l*Q_z + m
    if return_values:
        return pd.DataFrame({   'order'     :f90sussix.jklm.order_vec[:number_of_harmonics].astype(int),
                                'j'         :f90sussix.jklm.j_vec[:number_of_harmonics].astype(int),
                                'k'         :f90sussix.jklm.k_vec[:number_of_harmonics].astype(int),
                                'l'         :f90sussix.jklm.l_vec[:number_of_harmonics].astype(int),
                                'm'         :f90sussix.jklm.m_vec[:number_of_harmonics].astype(int),
                                'ax'        :kwargs['ax'][:number_of_harmonics],
                                'fx'        :np.angle(f90sussix.fcoe.zxpes[:number_of_harmonics]),
                                'ox'        :kwargs['ox'][:number_of_harmonics],
                                'ay'        :kwargs['ay'][:number_of_harmonics],
                                'fy'        :np.angle(f90sussix.fcoe.zypes[:number_of_harmonics]),
                                'oy'        :kwargs['oy'][:number_of_harmonics],
                                'as'        :kwargs['as'][:number_of_harmonics],
                                'fs'        :np.angle(f90sussix.fcoe.zspes[:number_of_harmonics]),
                                'os'        :kwargs['os'][:number_of_harmonics]})
        
