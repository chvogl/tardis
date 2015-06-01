from tardis.atomic.plugins import BaseAtomicDataType
import pandas as pd
from astropy import units as u
from collections import OrderedDict
import numpy as np
import h5py

class AtomicData(object):
    # all_data_types = BaseAtomicDataType.__subclasses__()

    @classmethod
    def from_hdf(cls, fname):
        return cls(pd.HDFStore(fname, mode='r'))

    def __init__(self, hdf_store, data_types=BaseAtomicDataType.__subclasses__()):
        self.data_types = [data_type()
                           for data_type in data_types]

        for plugin in self.data_types:
            try:
                plugin.load_hdf(hdf_store)
            except:
                print("Fail! at {0}".format(plugin.hdf_name))

        self.name2data_type = {item.hdf_name: item for item in self.data_types}



    def __getattr__(self, item):
        return self.name2data_type[item].value


class AtomicDataLeg(AtomicData):
    """
 attributes of `AtomData`   |  implementation status
 _____________________________________________________
 '_levels'                  |  implementation started (at the moment still problems with indexing in new `AtomicData`)
 'has_synpp_refs'           |  not implemented (preliminary implementation similar to has_macro_atom possible)
 'symbol2atomic_number'     |  implemented
 'has_ion_cx_data'          |  not implemented
 'ionization_data'          |  implemented
 'atomic_number2symbol'     |  implemented
 'macro_atom_references_all'|  not implemented
 'synpp_refs'               |  not implemented
 'data_sources'             |  ? not needed
 'macro_atom_data_all'      |  not implemented
 'version'                  |  not needed
 'zeta_data'                |  not implemented
 'has_zeta_data'            |  not implemented (preliminary implementation similar to has_macro_atom possible)
 'uuid1'                    |  not needed
 'has_macro_atom'           |  preliminary implementation
 'has_collision_data'       |  not implemented
 'md5'                      |  not needed
 '_lines'                   |  not implemented
 'atom_data'                |  implemented (however new 'AtomicData' has no masses in `atoms` at the moment)

  Functionality of the `prepare_atom_data` method of the old `AtomData` has to be implemented.
    """

    _atom_data = None
    _ionization_data = None
    _atomic_number2symbol = None
    _symbol2atomic_number = None
    __levels = None
    __lines = None
    _synpp_refs = None
    _macro_atom_data_all = None
    _zeta_data = None
    _has_macro_atom = None
    _has_collision_data = None
    _has_zeta_data = None
    _has_synpp_refs = None
    _has_ion_cx_data = None

    def __init__(self, hdf_store, data_types=BaseAtomicDataType.__subclasses__()):
        super(AtomicDataLeg, self).__init__(hdf_store, data_types)

    @property
    def atom_data(self):
        if self._atom_data is None:
            self._atom_data = self.atoms.copy()
            del self._atom_data['group']
            del self._atom_data['id']
            del self._atom_data['period']
            del self._atom_data['atomic_number']
            self._atom_data.set_index(self._atom_data.index + 1, inplace=True)
            cols = self._atom_data.columns.tolist()
            cols = [cols[i] for i in [1, 0]]
            self._atom_data = self._atom_data[cols]
            self._atom_data.index.rename('atomic_number', inplace=True)
        return self._atom_data

    @property
    def ionization_data(self):
        if self._ionization_data is None:
            self._ionization_data = self.ions.copy()
            self._ionization_data['ion_number'] += 1
            self._ionization_data['ionization_energy'] *= u.eV.to(u.erg)
            self._ionization_data.set_index(['atomic_number', 'ion_number'], inplace=True)
        return self._ionization_data

    @property
    def atomic_number2symbol(self):
        if self._atomic_number2symbol is None:
            self._atomic_number2symbol = OrderedDict(zip(np.array(self.atoms['atomic_number']),
                                                         np.array(self.atoms['symbol'])))
        return self._atomic_number2symbol

    @property
    def symbol2atomic_number(self):
        if self._symbol2atomic_number is None:
            self._symbol2atomic_number = OrderedDict(zip(np.array(self.atoms['symbol']),
                                                         np.array(self.atoms['atomic_number'])))
        return self._symbol2atomic_number

    @property
    def _levels(self):
        if self.__levels is None:
            self.__levels = self.levels.copy()
            cols = self.__levels.columns.tolist()
            cols = [cols[i] for i in [0, 3, 4, 1, 2]]
            self.__levels = self.__levels[cols]
            self.__levels['energy'] *= u.eV.to(u.erg)
        return self.__levels

    @property
    def _lines(self):
        if self.__lines is None:
            pass
        return self.__lines

    @property
    def md5(self):
        raise NotImplementedError

    @property
    def uuid1(self):
        raise NotImplementedError

    @property
    def version(self):
        raise NotImplementedError

    @property
    def has_macro_atom(self):
        """Has to be adjusted depending on the implementation of macro atom data in the new `AtomicData`.
         Here an `hdf-name` macro_atom is assumed for the `class MacroAtom(BaseAtomicDataType)` in plugins.py"""
        if self._has_macro_atom is None:
            if hasattr(self, 'macro_atom'):
                self._has_macro_atom = True
            else:
                self._has_macro_atom = False
        return self._has_macro_atom

    @atom_data.setter
    def atom_data(self, value):
        """
        This function returns the basic atomic data in the same way as the old atomic.py


        :param value:
        :return: numpy
        """
        # Do stuff with value
        self._atom_data = value


