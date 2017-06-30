"""
This is going to be the basic class which all classes inherit from.
"""
import numpy as np




class Base(object):

    def __init__(self, d_data):
        """ Init analysis object from given data dictionary. """
        self.data_id = d_data['data_id']
        self.meta = d_data['meta']
        self.spot = self.create_spot(gal_vir=d_data["gal_vir"],
                                     gal_pep=d_data["gal_pep"],
                                     data=d_data.get("data", None),
                                     data_std=d_data.get("data_std", None))

    @staticmethod
    def create_spot(gal_vir, gal_pep, data=None, data_std=None):
        """ Creates pandas DataFrame from input data.
        Subsequent analysis use the spot DataFrame.

        :return:
        """
        vir_cor = gal_vir.pivot(index="Row", columns="Column", values="Name")
        pep_cor = gal_pep.pivot(index="Row", columns="Column", values="Name")

        # merge complete spotinformation
        vir_cor_unstacked = vir_cor.unstack()
        spot = pep_cor.unstack()
        spot = spot.reset_index()
        spot = spot.rename(columns={0: "Peptide"})
        spot["Reference"] = spot["Peptide"].str.contains("DYE")
        spot["Virus"] = vir_cor_unstacked.values
        spot["Buffer"] = spot["Peptide"].str.match("Buffer")
        spot["No_Peptide"] = spot["Peptide"].str.match("NO")



        if data is None:
            spot["Intensity"] = np.NaN
        else:
            spot["Intensity"] = data.unstack().values

        if data_std is None:
            spot["Std"] = np.NaN
        else:
            spot["Std"] = data_std.unstack().values


        # how often this spot measured
        spot["Replica"] = 0
        for virus_unique in spot["Virus"].unique():
            for peptide_unique in spot["Peptide"].unique():
                replica = 0
                for index in spot.index:
                    if spot["Virus"][index] == virus_unique and spot["Peptide"][index] == peptide_unique:
                        spot.set_value(index, "Replica", replica)
                        replica += 1

        return spot