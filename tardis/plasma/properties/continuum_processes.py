import logging

import numpy as np
import pandas as pd

from numba import prange, njit
from astropy import constants as const

from tardis.plasma.exceptions import PlasmaException
from tardis.plasma.properties.base import (ProcessingPlasmaProperty,
                                           Input)
from tardis.plasma.properties.j_blues import JBluesDiluteBlackBody

__all__ = ['SpontRecombRateCoeff', 'StimRecombRateCoeff', 'PhotoIonRateCoeff',
           'PhotoIonEstimatorsNormFactor', 'PhotoIonRateCoeffEstimator',
           'StimRecombRateCoeffEstimator', 'CorrPhotoIonRateCoeff',
           'BfHeatingRateCoeffEstimator', 'SpontRecombCoolingRateCoeff',
           'BaseRecombTransProbs', 'BasePhotoIonTransProbs',
           'CollDeexcRateCoeff', 'CollExcRateCoeff', 'BaseCollisionTransProbs']

logger = logging.getLogger(__name__)

njit_dict = {'fastmath': False, 'parallel': False}


@njit(**njit_dict)
def integrate_array_by_blocks(f, x, block_references):
    """
    Integrates a function f defined at locations x over blocks
    given in block_references.

    Parameters
    ----------
    f : Two-dimensional Numpy Array, dtype float
    x : One-dimensional Numpy Array, dtype float
    block_references : One-dimensional Numpy Array, dtype int

    Returns
    -------
    integrated : Two-dimensional Numpy Array, dtype float

    """
    integrated = np.zeros((len(block_references) - 1, f.shape[1]))
    for i in prange(f.shape[1]):  # columns
        for j in prange(len(integrated)):  # rows
            start = block_references[j]
            stop = block_references[j + 1]
            integrated[j, i] = np.trapz(f[start:stop, i], x[start:stop])
    return integrated


def get_ion_multi_index(multi_index_full, next_higher=True):
    """
    Calculates the corresponding ion MultiIndex for a level MultiIndex.

    Parameters
    ----------
    multi_index_full : Pandas MultiIndex (atomic_number, ion_number,
                                          level_number)
    next_higher : bool
        If true use ion number of next higher ion, else use ion_number from
        multi_index_full.

    Returns
    -------
    multi_index : Pandas MultiIndex (atomic_number, ion_number)

    """
    atomic_number = multi_index_full.get_level_values(0)
    ion_number = multi_index_full.get_level_values(1)
    if next_higher is True:
        ion_number += 1
    return pd.MultiIndex.from_arrays([atomic_number, ion_number])


def get_ground_state_multi_index(multi_index_full):
    """
    Calculates the level MultiIndex for the ground state of the next higher
    ion.

    Parameters
    ----------
    multi_index_full : Pandas MultiIndex (atomic_number, ion_number,
                                          level_number)

    Returns
    -------
    multi_index : Pandas MultiIndex (atomic_number, ion_number)

    """
    atomic_number = multi_index_full.get_level_values(0)
    ion_number = multi_index_full.get_level_values(1) + 1
    level_number = np.zeros_like(ion_number)
    return pd.MultiIndex.from_arrays([atomic_number, ion_number, level_number])


class IndexSetterMixin(object):
    @staticmethod
    def set_index(p, photo_ion_idx, transition_type=0, reverse=True):
        idx = photo_ion_idx.loc[p.index]
        transition_type = transition_type * np.ones_like(
            idx.destination_level_idx
        )
        transition_type = pd.Series(transition_type, name='transition_type')
        idx_arrays = [idx.source_level_idx, idx.destination_level_idx]
        if reverse:
            idx_arrays = idx_arrays[::-1]
        idx_arrays.append(transition_type)
        index = pd.MultiIndex.from_arrays(idx_arrays)
        if reverse:
            index.names = index.names[:-1][::-1] + [index.names[-1]]
        p = p.set_index(index, drop=True)
        return p


class SpontRecombRateCoeff(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    alpha_sp : Pandas DataFrame, dtype float
               The rate coefficient for spontaneous recombination.
    """
    outputs = ('alpha_sp',)
    latex_name = ('\\alpha^{\\textrm{sp}}',)

    def calculate(self, photo_ion_cross_sections, t_electrons,
                  photo_ion_block_references, photo_ion_index, phi_ik):
        x_sect = photo_ion_cross_sections['x_sect'].values
        nu = photo_ion_cross_sections['nu'].values

        alpha_sp = (8 * np.pi * x_sect * nu ** 2 / (const.c.cgs.value) ** 2)
        alpha_sp = alpha_sp[:, np.newaxis]
        boltzmann_factor = np.exp(-nu[np.newaxis].T / t_electrons *
                                  (const.h.cgs.value / const.k_B.cgs.value))
        alpha_sp = alpha_sp * boltzmann_factor
        alpha_sp = integrate_array_by_blocks(alpha_sp, nu,
                                             photo_ion_block_references)
        alpha_sp = pd.DataFrame(alpha_sp, index=photo_ion_index)
        return alpha_sp * phi_ik.loc[alpha_sp.index]


class SpontRecombCoolingRateCoeff(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    c_fb_sp : Pandas DataFrame, dtype float
              The rate coefficient for cooling by
              spontaneous recombination.
    """
    outputs = ('c_fb_sp',)
    latex_name = ('ca^{\\textrm{sp}}_{\\textrm{fb}}',)

    def calculate(self, photo_ion_cross_sections, t_electrons,
                  photo_ion_block_references, photo_ion_index, phi_ik, nu_i):
        x_sect = photo_ion_cross_sections['x_sect'].values
        nu = photo_ion_cross_sections['nu'].values
        factor = (1 - nu_i / photo_ion_cross_sections['nu']).values
        alpha_sp = (8 * np.pi * x_sect * factor * nu ** 3 /
                    (const.c.cgs.value) ** 2)
        alpha_sp = alpha_sp[:, np.newaxis]
        boltzmann_factor = np.exp(-nu[np.newaxis].T / t_electrons *
                                  (const.h.cgs.value / const.k_B.cgs.value))
        alpha_sp = alpha_sp * boltzmann_factor
        alpha_sp = integrate_array_by_blocks(alpha_sp, nu,
                                             photo_ion_block_references)
        alpha_sp = pd.DataFrame(alpha_sp, index=photo_ion_index)
        return alpha_sp * phi_ik.loc[alpha_sp.index]


class PhotoIonRateCoeff(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    gamma : Pandas DataFrame, dtype float
            The rate coefficient for radiative ionization.
    """
    outputs = ('gamma',)
    latex_name = ('\\gamma',)

    def calculate(self, photo_ion_cross_sections, gamma_estimator,
                  photo_ion_norm_factor, photo_ion_block_references,
                  photo_ion_index, t_rad, w):
        # Used for initialization
        if gamma_estimator is None:
            gamma = self.calculate_from_dilute_bb(photo_ion_cross_sections,
                                                  photo_ion_block_references,
                                                  photo_ion_index, t_rad, w)
        else:
            gamma = gamma_estimator * photo_ion_norm_factor
        return gamma

    @staticmethod
    def calculate_from_dilute_bb(photo_ion_cross_sections,
                                 photo_ion_block_references,
                                 photo_ion_index, t_rad, w):
        nu = photo_ion_cross_sections['nu']
        x_sect = photo_ion_cross_sections['x_sect']
        j_nus = JBluesDiluteBlackBody.calculate(
            photo_ion_cross_sections, nu, t_rad, w
        )
        gamma = j_nus.multiply(
            4. * np.pi * x_sect / nu / const.h.cgs.value, axis=0
        )
        gamma = integrate_array_by_blocks(gamma.values, nu.values,
                                          photo_ion_block_references)
        gamma = pd.DataFrame(gamma, index=photo_ion_index)
        return gamma


class StimRecombRateCoeff(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    alpha_stim : Pandas DataFrame, dtype float
                 The rate coefficient for stimulated recombination.
    """
    outputs = ('alpha_stim',)
    latex_name = ('\\alpha^{\\textrm{stim}}',)

    def calculate(self, photo_ion_cross_sections, alpha_stim_estimator,
                  photo_ion_norm_factor, photo_ion_block_references,
                  photo_ion_index, t_rad, w, phi_ik, t_electrons):
        # Used for initialization
        if alpha_stim_estimator is None:
            alpha_stim = self.calculate_from_dilute_bb(
                photo_ion_cross_sections, photo_ion_block_references,
                photo_ion_index, t_rad, w, t_electrons
            )
            alpha_stim *= phi_ik.loc[alpha_stim.index]
        else:
            alpha_stim = alpha_stim_estimator * photo_ion_norm_factor
        return alpha_stim

    @staticmethod
    def calculate_from_dilute_bb(photo_ion_cross_sections,
                                 photo_ion_block_references,
                                 photo_ion_index, t_rad, w, t_electrons):
        nu = photo_ion_cross_sections['nu']
        x_sect = photo_ion_cross_sections['x_sect']
        boltzmann_factor = np.exp(-nu.values[np.newaxis].T / t_electrons *
                                  (const.h.cgs.value / const.k_B.cgs.value))
        j_nus = JBluesDiluteBlackBody.calculate(
            photo_ion_cross_sections, nu, t_rad, w
        )
        j_nus *= boltzmann_factor
        alpha_stim = j_nus.multiply(
            4. * np.pi * x_sect / nu / const.h.cgs.value, axis=0
        )
        alpha_stim = integrate_array_by_blocks(alpha_stim.values, nu.values,
                                               photo_ion_block_references)
        alpha_stim = pd.DataFrame(alpha_stim, index=photo_ion_index)
        return alpha_stim


class BaseRecombTransProbs(ProcessingPlasmaProperty, IndexSetterMixin):
    """
    Attributes
    ----------
    p_recomb : Pandas DataFrame, dtype float
               The unnormalized transition probabilities for
               spontaneous recombination.
    """
    outputs = ('p_recomb', )
    latex_name = ('p^{\\textrm{recomb}}', '')

    def calculate(self, alpha_sp, nu_i, energy_i, photo_ion_idx):
        p_recomb_deac = alpha_sp.multiply(nu_i, axis=0) * const.h.cgs.value
        p_recomb_deac = self.set_index(p_recomb_deac, photo_ion_idx,
                                       transition_type=-1)

        p_recomb_internal = alpha_sp.multiply(energy_i, axis=0)
        p_recomb_internal = self.set_index(p_recomb_internal, photo_ion_idx,
                                           transition_type=0)
        p_recomb = pd.concat([p_recomb_deac, p_recomb_internal])
        return p_recomb


class BasePhotoIonTransProbs(ProcessingPlasmaProperty, IndexSetterMixin):
    """
    Attributes
    ----------
    p_photo_ion : Pandas DataFrame, dtype float
                  The unnormalized transition probabilities for
                  radiative ionization.
    """
    outputs = ('p_photo_ion', )
    latex_name = ('p^{\\textrm{photo_ion}}', )

    def calculate(self, gamma_corr, nu_i, photo_ion_idx):
        p_photo_ion = gamma_corr.multiply(nu_i, axis=0) * const.h.cgs.value
        p_photo_ion = self.set_index(p_photo_ion, photo_ion_idx, reverse=False)
        return p_photo_ion


class CorrPhotoIonRateCoeff(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    gamma : Pandas DataFrame, dtype float
            The rate coefficient for radiative ionization corrected for
            stimulated recombination.
    """
    outputs = ('gamma_corr',)
    latex_name = ('\\gamma_\\mathrm{corr}',)

    def calculate(self, gamma, alpha_stim, electron_densities,
                  ion_number_density, level_number_density):
        n_k_index = get_ion_multi_index(alpha_stim.index)
        n_k = ion_number_density.loc[n_k_index].values
        n_i = level_number_density.loc[alpha_stim.index].values
        gamma_corr = gamma - alpha_stim * n_k * electron_densities / n_i
        num_neg_elements = (gamma_corr < 0).sum().sum()
        if num_neg_elements:
            raise PlasmaException(
                "Negative values in CorrPhotoIonRateCoeff."
            )
        return gamma_corr


class PhotoIonEstimatorsNormFactor(ProcessingPlasmaProperty):
    outputs = ('photo_ion_norm_factor',)
    latex = ('\\frac{1}}}{'
             'time_\\textrm{simulation} volume h}')

    @staticmethod
    def calculate(time_simulation, volume):
        return (time_simulation * volume * const.h.cgs.value)**-1


class PhotoIonRateCoeffEstimator(Input):
    """
    Attributes
    ----------
    gamma_estimator : Pandas DataFrame, dtype float
                      Unnormalized MC estimator for the rate coefficient
                      for radiative ionization.
    """
    outputs = ('gamma_estimator',)


class StimRecombRateCoeffEstimator(Input):
    """
    Attributes
    ----------
    alpha_stim_estimator : Pandas DataFrame, dtype float
                           Unnormalized MC estimator for the rate coefficient
                           for stimulated recombination.
    """
    outputs = ('alpha_stim_estimator',)


class BfHeatingRateCoeffEstimator(Input):
    """
    Attributes
    ----------
    bf_heating_coeff_estimator : Pandas DataFrame, dtype float
                                 Unnormalized MC estimator for the rate
                                 coefficient for bound-free heating.
    """
    outputs = ('bf_heating_coeff_estimator',)


class CollExcRateCoeff(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    coll_exc_coeff : Pandas DataFrame, dtype float
        Rate coefficient for collisional excitation.
    """
    outputs = ('coll_exc_coeff',)
    latex_name = ('c_{lu}',)

    def calculate(self, yg_interp, yg_index, t_electrons, delta_E_yg):
        yg = yg_interp(t_electrons)
        k_B = const.k_B.cgs.value
        boltzmann_factor = np.exp(
            - delta_E_yg.values[np.newaxis].T / (t_electrons * k_B)
        )
        q_ij = 8.629e-6 / np.sqrt(t_electrons) * yg * boltzmann_factor
        return pd.DataFrame(q_ij, index=yg_index)


class CollDeexcRateCoeff(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    coll_deexc_coeff : Pandas DataFrame, dtype float
        Rate coefficient for collisional deexcitation.
    """
    outputs = ('coll_deexc_coeff',)
    latex_name = ('c_{ul}',)

    def calculate(self, thermal_lte_level_boltzmann_factor, coll_exc_coeff):
        ll_index = coll_exc_coeff.index.droplevel('level_number_upper')
        lu_index = coll_exc_coeff.index.droplevel('level_number_lower')

        n_lower_prop = thermal_lte_level_boltzmann_factor.loc[ll_index].values
        n_upper_prop = thermal_lte_level_boltzmann_factor.loc[lu_index].values

        coll_deexc_coeff = coll_exc_coeff * n_lower_prop / n_upper_prop
        return coll_deexc_coeff


class BaseCollisionTransProbs(ProcessingPlasmaProperty, IndexSetterMixin):
    """
    Attributes
    ----------
    p_coll : Pandas DataFrame, dtype float
        The unnormalized transition probabilities for
        collisional excitation.
    """
    outputs = ('p_coll', )
    latex_name = ('p^{\\textrm{coll}}', '')

    def calculate(self, coll_exc_coeff, coll_deexc_coeff, yg_idx,
                  electron_densities, delta_E_yg, atomic_data,
                  level_number_density):
        p_deexc_deac = (coll_deexc_coeff * electron_densities).multiply(
            delta_E_yg.values, axis=0)
        p_deexc_deac = self.set_index(p_deexc_deac, yg_idx)
        p_deexc_deac = p_deexc_deac.groupby(level=[0]).sum()
        index_dd = pd.MultiIndex.from_product(
            [p_deexc_deac.index.values, ['k'], [0]],
            names=list(yg_idx.columns) + ['transition_type']
        )
        p_deexc_deac = p_deexc_deac.set_index(index_dd)

        ll_index = coll_deexc_coeff.index.droplevel('level_number_upper')
        energy_lower = atomic_data.levels.energy.loc[ll_index]
        p_deexc_internal = (coll_deexc_coeff * electron_densities).multiply(
            energy_lower.values, axis=0)
        p_deexc_internal = self.set_index(p_deexc_internal, yg_idx,
                                          transition_type=0, reverse=True)

        p_exc_internal = (coll_exc_coeff * electron_densities).multiply(
            energy_lower.values, axis=0)
        p_exc_internal = self.set_index(p_exc_internal, yg_idx,
                                        transition_type=0, reverse=False)
        p_exc_cool = (coll_exc_coeff * electron_densities).multiply(
            delta_E_yg.values, axis=0)
        p_exc_cool = p_exc_cool * level_number_density.loc[ll_index].values
        p_exc_cool = self.set_index(p_exc_cool, yg_idx, reverse=False)
        p_exc_cool = p_exc_cool.groupby(level='destination_level_idx').sum()
        exc_cool_index = pd.MultiIndex.from_product(
            [['k'], p_exc_cool.index.values, [0]],
            names=list(yg_idx.columns) + ['transition_type']
        )
        p_exc_cool = p_exc_cool.set_index(exc_cool_index)
        p_coll = pd.concat(
            [p_deexc_deac, p_deexc_internal, p_exc_internal, p_exc_cool]
        )
        return p_coll
