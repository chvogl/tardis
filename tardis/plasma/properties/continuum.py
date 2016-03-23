import logging

import numpy as np
import pandas as pd
from scipy.integrate import trapz
from astropy import units as u, constants as const

from tardis.plasma.properties.base import ProcessingPlasmaProperty

logger = logging.getLogger(__name__)

__all__ = ['SpontRecombRateCoeff']

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