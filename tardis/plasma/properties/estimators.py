from tardis.plasma.properties.base import (Input, ArrayInput, DataFrameInput)

__all__ = ['PhotoIonRateEstimator', 'StimRecombRateEstimator', 'PhotoIonRateStatistics']

class PhotoIonRateEstimator(DataFrameInput):
    """
    Attributes
    ----------
    photo_ion_estimator : Pandas DataFrame
                          Monte Carlo estimator for the photoionization rate coefficient.
    """
    outputs = ('photo_ion_estimator',)
    latex_name = ('\\gamma_{\\textrm{estimator}}',)

class PhotoIonRateStatistics(DataFrameInput):
    """
    Attributes
    ----------
    photo_ion_statistics : Pandas DataFrame
                           Number of updates of the photoionization rate estimator.
    """
    outputs = ('photo_ion_statistics',)
    latex_name = ('\\N_{\\textrm{bf_estimator}}',)

class StimRecombRateEstimator(DataFrameInput):
    """
    Attributes
    ----------
    stim_recomb_estimator : Pandas DataFrame
                            Monte Carlo estimator for the stimulated recombination rate coefficient.
    """
    outputs = ('stim_recomb_estimator',)
    latex_name = ('\\alpha_{\\textrm{estimator}}^{\\textrm{st}}',)