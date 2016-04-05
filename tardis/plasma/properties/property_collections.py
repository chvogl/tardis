from tardis.plasma.properties import *
from tardis.plasma.properties.estimators import *
from tardis.plasma.properties.continuum import *

class PlasmaPropertyCollection(list):
    pass

basic_inputs = PlasmaPropertyCollection([TRadiative, Abundance, Density,
    TimeExplosion, AtomicData, JBlues, DilutionFactor, LinkTRadTElectron])
basic_properties = PlasmaPropertyCollection([BetaRadiation,
    Levels, Lines, AtomicMass, PartitionFunction,
    GElectron, IonizationData, NumberDensity, LinesLowerLevelIndex,
    LinesUpperLevelIndex, TauSobolev, LevelNumberDensity,
    StimulatedEmissionFactor, SelectedAtoms, ElectronTemperature])
non_nlte_ionzation_properties = PlasmaPropertyCollection([IonNumberDensity])
lte_ionization_properties = PlasmaPropertyCollection([PhiSahaLTE])
lte_excitation_properties = PlasmaPropertyCollection([LevelBoltzmannFactorLTE])
macro_atom_properties = PlasmaPropertyCollection([BetaSobolev,
    TransitionProbabilities])
nebular_ionization_properties = PlasmaPropertyCollection([PhiSahaNebular,
    ZetaData, BetaElectron, RadiationFieldCorrection])
continuum_inputs = PlasmaPropertyCollection([StimRecombRateEstimator, PhotoIonRateEstimator,
    PhotoIonRateStatistics, ContinuumData])
continuum_lte_properties = PlasmaPropertyCollection([
    PhiSahaLTECont, LTEIonNumberDensity, LevelBoltzmannFactorLTECont, LTEPartitionFunction, LTELevelNumberDensity,
    PhiSahaElectrons, GElectronTe, LevelBoltzmannFactorLTETe, LTEPartitionFunctionTe, PhiLucy])
continuum_interaction_properties = PlasmaPropertyCollection([
    SpontRecombRateCoeff, StimRecombRateCoeff,PhotoIonRateCoeff])
nlte_ionizaton_properties = PlasmaPropertyCollection([NLTEIonNumberDensity])
dilute_lte_excitation_properties = PlasmaPropertyCollection([
    LevelBoltzmannFactorDiluteLTE])
non_nlte_properties = PlasmaPropertyCollection([LevelBoltzmannFactorNoNLTE])
nlte_properties = PlasmaPropertyCollection([
    LevelBoltzmannFactorNLTE, NLTEData, LTEJBlues, PreviousElectronDensities,
    PreviousBetaSobolev, BetaSobolev])
helium_nlte_properties = PlasmaPropertyCollection([HeliumNLTE,
    RadiationFieldCorrection, ZetaData,
    BetaElectron])
helium_numerical_nlte_properties = PlasmaPropertyCollection([
    HeliumNumericalNLTE])