from tardis.atomic.plugins import BaseAtomicDataType
import pandas as pd

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
    _atom_data = None

    def __init__(self, hdf_store, data_types):
        super_object = super(AtomicDataLeg, self)
        super_object.__init__(hdf_store, data_types)


    @property
    def atom_data(self):
        return  # The new atom data with old shape


    @atom_data.setter
    def atom_data(self, value):
        """
        This function returns the basic atomic data in the same way as the old atomic.py


        :param value:
        :return: numpy
        """
        # Do stuff with value
        self._atom_data = value


