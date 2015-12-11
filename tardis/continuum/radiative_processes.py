import logging

import numpy as np
from astropy import constants as const
import pandas as pd
from scipy.integrate import simps

from tardis.util import intensity_black_body
from tardis.continuum.base import PhysicalContinuumProcess, BoundFreeEnergyMixIn
from tardis.continuum.constants import continuum_constants as cconst


logger = logging.getLogger(__name__)


class RadiativeIonization(PhysicalContinuumProcess, BoundFreeEnergyMixIn):
    name = 'radiative_ionization'
    cooling = False
    macro_atom_transitions = 'continuum'

    def __init__(self, input_data):
        super(RadiativeIonization, self).__init__(input_data)

    def _calculate_rate_coefficient(self):
        # Corrected photoionization coefficient
        j_nus = self._calculate_j_nus()
        stimulated_emission_correction = self._calculate_stimulated_emission_correction()
        corrected_photoion_coeff = j_nus.multiply(4. * np.pi * self.photoionization_data['x_sect'] /
                                                  self.photoionization_data['nu'] / const.h.cgs.value, axis=0)
        corrected_photoion_coeff = corrected_photoion_coeff.multiply(stimulated_emission_correction)
        corrected_photoion_coeff.insert(0, 'nu', self.photoionization_data['nu'])
        corrected_photoion_coeff = corrected_photoion_coeff.groupby(level=[0, 1, 2])
        tmp = {}
        for i in range(len(self.t_rads)):
            tmp[i] = corrected_photoion_coeff.apply(lambda sub: simps(sub[i], sub['nu'], even='first'))
        corrected_photoion_coeff = pd.DataFrame(tmp)
        return corrected_photoion_coeff

    def _calculate_j_nus(self):
        nus = self.photoionization_data['nu'].values
        j_nus = self.ws * intensity_black_body(nus[np.newaxis].T, self.t_rads)
        return pd.DataFrame(j_nus, index=self.photoionization_data.index, columns=np.arange(len(self.t_rads)))

    def _calculate_boltzmann_factor(self, nu):
        u0s = self._calculate_u0s(nu)
        return np.exp(-u0s)

    def _calculate_stimulated_emission_correction(self):
        nu = self.photoionization_data['nu'].values
        boltzmann_factor = self._calculate_boltzmann_factor(nu)
        lte_nonlte_level_pop_ratio = self._get_lte_nonlte_level_pop_ratio(self.photoionization_data.index)
        correction_factor = (1. - lte_nonlte_level_pop_ratio * boltzmann_factor)
        return correction_factor

    def _get_lte_nonlte_level_pop_ratio(self, index):
        level_pop_lte = self._get_lte_level_pop(index)
        level_pop = self._get_level_pop(index)
        continuum_pop = self._get_ion_number_density(index)
        continuum_pop_lte = self._get_lte_ion_number_density(index)
        ratio = (continuum_pop / continuum_pop_lte) * (level_pop_lte / level_pop)
        return ratio

    @property
    def level_lower_energy(self):
        return self._get_level_energy(self.rate_coefficient.index)


class RadiativeRecombination(PhysicalContinuumProcess, BoundFreeEnergyMixIn):
    name = 'radiative_recombination'

    def __init__(self, input_data):
        super(RadiativeRecombination, self).__init__(input_data)

    def _calculate_cooling_rate(self):
        sp_recombination_coeff_E = self._calculate_rate_coefficient(modified=True)
        fb_cooling_rate = (sp_recombination_coeff_E - self.rate_coefficient)
        fb_cooling_rate = fb_cooling_rate.multiply(const.h.cgs.value * self.nu_i, axis=0)
        fb_cooling_rate = fb_cooling_rate.multiply(self.electron_densities, axis=1)
        ion_number_density = self._get_ion_number_density(fb_cooling_rate.index)
        fb_cooling_rate = fb_cooling_rate.multiply(ion_number_density)
        continuum_edge_idx = self._get_continuum_edge_idx(fb_cooling_rate.index)
        fb_cooling_rate.set_index(continuum_edge_idx, inplace=True)
        return fb_cooling_rate

    def _calculate_rate_coefficient(self, modified=False):
        # Spontaneous recombination coefficient
        if modified == False:
            recomb_coeff = (8 * np.pi * self.photoionization_data['x_sect']
                            * (self.photoionization_data['nu']) ** 2 / (const.c.cgs.value) ** 2).values
        else:
            recomb_coeff = (8 * np.pi * self.photoionization_data['x_sect']
                            * (self.photoionization_data['nu']) ** 3 / (const.c.cgs.value) ** 2).values

        recomb_coeff = recomb_coeff[:, np.newaxis]
        boltzmann_factor = np.exp(-self.photoionization_data.nu.values[np.newaxis].T / \
                                  self.t_rads * (const.h.cgs.value / const.k_B.cgs.value))
        recomb_coeff = pd.DataFrame(boltzmann_factor * recomb_coeff, index=self.photoionization_data.index)
        recomb_coeff = recomb_coeff.divide(self.electron_densities, axis=1)
        recomb_coeff.insert(0, 'nu', self.photoionization_data['nu'])
        recomb_coeff = recomb_coeff.groupby(level=[0, 1, 2])
        tmp = {}
        for i in range(self.no_of_shells):
            tmp[i] = recomb_coeff.apply(lambda sub: simps(sub[i], sub['nu'], even='first'))
            if modified == True:
                tmp[i] /= self.nu_i
        recomb_coeff = pd.DataFrame(tmp)
        recomb_coeff = recomb_coeff.multiply(self._get_lte_level_pop(recomb_coeff.index))
        ion_number_density = self._get_ion_number_density(recomb_coeff.index, dtype='dataframe')
        recomb_coeff = recomb_coeff.divide(ion_number_density.values)
        return recomb_coeff


class RadiativeExcitation(PhysicalContinuumProcess):
    name = 'radiative_excitation'
    cooling = False
    macro_atom_transitions = 'up'

    @property
    def internal_jump_probabilities(self):
        return self.input.radiative_transition_probabilities_prep.loc[self.transition_up_filter] * cconst.c_einstein


class RadiativeDeexcitation(PhysicalContinuumProcess):
    name = 'radiative_deexcitation'
    cooling = False
    macro_atom_transitions = 'down'

    @property
    def internal_jump_probabilities(self):
        return self.input.radiative_transition_probabilities_prep.loc[self.transition_down_filter] * cconst.c_einstein

    @property
    def deactivation_probabilities(self):
        filter = self.transition_deactivation_filter
        deactivation_probabilities = self.input.radiative_transition_probabilities_prep.loc[filter] * cconst.c_einstein
        deactivation_probabilities.insert(0, 'lines_idx', self.macro_atom_data.loc[filter, 'lines_idx'].values)
        return deactivation_probabilities


class FreeFree(PhysicalContinuumProcess):
    name = 'free_free'

    def _calculate_cooling_rate(self, **kwargs):
        # TODO: value for Gaunt factor (Lucy: = 1; Osterbrock recommendation for nebular conditions: = 1.3 )
        factor = self.ion_number_density.mul(np.square(self.input.ion_charges), axis=0).sum().values
        cooling_rate = cconst.C0_ff * self.electron_densities * np.sqrt(self.t_electrons) * factor
        return cooling_rate