import logging

import numpy as np
import pandas as pd
from scipy.integrate import trapz
from astropy import units as u, constants as const

from tardis.plasma.properties.base import ProcessingPlasmaProperty

logger = logging.getLogger(__name__)

__all__ = ['SpontRecombRateCoeff', 'PhotoIonRateCoeff', 'StimRecombRateCoeff']

def get_estimator_index(photo_ion_index_sorted):
    index = pd.MultiIndex.from_tuples(photo_ion_index_sorted)
    index.names = [u'atomic_number', u'ion_number', u'level_number']
    return index

def calculate_rate_coefficient_from_estimator(estimator, photo_ion_index_sorted):
    no_of_shells = estimator.shape[1]
    index = get_estimator_index(photo_ion_index_sorted)
    rate_coeff = pd.DataFrame(estimator, index=index, columns=np.arange(no_of_shells))
    return rate_coeff

class SpontRecombRateCoeff(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    alpha_sp : Pandas DataFrame, dtype float
               The rate coefficient for spontaneous recombination.
    """
    outputs = ('alpha_sp',)
    latex_name = ('\\alpha^{\\textrm{sp}}',)
    # TODO
    latex_formula = ('',)

    def calculate(self, photo_ion_cross_sections, t_electrons, phi_lucy):
        x_sect = photo_ion_cross_sections['x_sect']
        nu = photo_ion_cross_sections['nu']

        alpha_sp = (8 * np.pi * x_sect * nu ** 2 / (const.c.cgs.value) ** 2).values
        alpha_sp = alpha_sp[:, np.newaxis]
        boltzmann_factor = np.exp(-nu.values[np.newaxis].T / t_electrons * (const.h.cgs.value / const.k_B.cgs.value))
        recomb_coeff = pd.DataFrame(boltzmann_factor * alpha_sp, index=nu.index)
        recomb_coeff.insert(0, 'nu', nu)
        recomb_coeff = recomb_coeff.groupby(level=[0, 1, 2])
        tmp = {}
        for i in range(len(t_electrons)):
            tmp[i] = recomb_coeff.apply(lambda sub: trapz(sub[i], sub['nu']))
        alpha_sp = pd.DataFrame(tmp)
        return alpha_sp.multiply(phi_lucy)

class PhotoIonRateCoeff(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    gamma : Pandas DataFrame, dtype float
               The rate coefficient for radiative ionization.
    """
    outputs = ('gamma',)
    latex_name = ('\\gamma',)
    # TODO
    latex_formula = ('',)

    def calculate(self, photo_ion_estimator, photo_ion_index_sorted):
        if photo_ion_estimator is not None:
            gamma = calculate_rate_coefficient_from_estimator(photo_ion_estimator, photo_ion_index_sorted)
        else:
            gamma = None
        return gamma

class StimRecombRateCoeff(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    gamma : Pandas DataFrame, dtype float
               The rate coefficient for radiative ionization.
    """
    outputs = ('alpha_stim',)
    latex_name = ('\\alpha^{\\textrm{stim}}',)
    # TODO
    latex_formula = ('',)

    def calculate(self, stim_recomb_estimator, photo_ion_index_sorted, phi_lucy):
        if stim_recomb_estimator is not None:
            alpha_stim = calculate_rate_coefficient_from_estimator(stim_recomb_estimator, photo_ion_index_sorted)
            alpha_stim = alpha_stim.multiply(phi_lucy)
        else:
            alpha_stim = None
        return alpha_stim