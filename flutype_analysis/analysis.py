"""
Analysis script.

Information for a given microtiter or microarray consists of the following files
with the following naming schema:
        {id}_pep.gal : peptide gal file (which peptide at which location)
        {id}_vir.gal : virus gal file (which virus at which location)
        {id}.meta : additional information
        {id}.csv : numerical data

The numerical data in the csv is either direct data from the fluorescence reader,
or spot intensities quantified from images.
"""
from __future__ import print_function, absolute_import
import os
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from sklearn import decomposition

from flutype_analysis import correlation
from flutype_analysis import utils


PATTERN_PEP_GAL = "{}_pep.gal"
PATTERN_VIR_GAL = "{}_vir.gal"
PATTERN_META = "{}.meta"
PATTERN_DATA = "{}.csv"
PATTERN_STD = "{}_std.csv"


def load_data(data_id, directory='data/flutype_test'):
    """ Loads datasets for given data_id.
    See naming schema above.

    :param data_id: id of the Data
    :return: dictionary of information
    """
    print("-" * 80)
    print("Loading data corresponding to data_id: <{}> in dir <{}>".format(data_id, directory))
    print("-" * 80)

    f_gal_vir = os.path.join(directory, PATTERN_VIR_GAL.format(data_id))
    gal_vir = pd.read_csv(f_gal_vir, sep='\t', index_col="ID")
    print("{}:{}".format("Virus .gal", f_gal_vir))

    f_gal_pep = os.path.join(directory, PATTERN_PEP_GAL.format(data_id))
    gal_pep = pd.read_csv(f_gal_pep, sep='\t', index_col="ID")
    print("Peptide .gal :{}".format(f_gal_pep))

    f_meta = os.path.join(directory, PATTERN_META.format(data_id))
    meta = pd.read_csv(f_meta, sep='\t')
    print("Meta  :{}".format(f_meta))

    f_data = os.path.join(directory, PATTERN_DATA.format(data_id))
    data = pd.read_csv(f_data, sep='\t', index_col=0)
    print("Spot intensity file  :{}".format(f_data))

    f_std = os.path.join(directory, PATTERN_STD.format(data_id))

    d = {'data_id': data_id,
         'data': data,
         'gal_pep': gal_pep,
         'gal_vir': gal_vir,
         'meta': meta}

    if os.path.exists(f_std):
        d['data_std'] = pd.read_csv(f_std, sep='\t', index_col=0)
        print("Intensity standard deviation:{}".format(f_std))
    else:
        print("Spot intensities for the data ID ({}) are not averaged but primary values".format(data_id))

    return d


class Analysis(object):

    def __init__(self, d_data):
        """ Init analysis object from given data dictionary. """
        self.data_id = d_data['data_id']
        self.meta = d_data['meta']

        self.spot = self.create_spot(gal_vir=d_data["gal_vir"],
                                     gal_pep=d_data["gal_pep"],
                                     data=d_data["data"],
                                     data_std=d_data.get("data_std", None))

    @staticmethod
    def create_spot(gal_vir, gal_pep, data, data_std=None):
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
        spot["Referenz"] = spot["Peptide"].str.contains("Leuchte")
        spot["Virus"] = vir_cor_unstacked.values
        spot["Intensity"] = data.unstack().values

        if data_std:
            spot["Std"] = data_std.unstack().values
        else:
            spot["Std"] = np.NaN

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

    def heatmap(self, descript=True, heatmap=True, **kwargs):
        """ Creates heatmap

        :param descript: boolean flag to show peptide information
        :param heatmap: boolean flag to show heatmap
        :param kwargs:
        :return:
        """

        fig, ax = plt.subplots(**kwargs)
        # imshow portion
        if descript:
            alpha = 0.5
        else:
            alpha = 1.0
        if heatmap:
            ax.imshow(self.spot.pivot(index='Column', columns='Row', values='Intensity'), interpolation='nearest',
                    cmap="hot", alpha=alpha)
        # plt.pcolor(Spot["Intensity"].unstack())
        # text portion
        x, y = np.meshgrid(self.spot["Row"].unique() - 1, self.spot["Column"].unique() - 1)
        for index in self.spot.index:
            # print(y_val)
            if descript:
                ax.text(self.spot["Row"].loc[index] - 0.90, self.spot["Column"].loc[index] - 0.90, self.spot["Peptide"].loc[index],
                        va='center', ha='center', fontsize=8, weight="bold", rotation=45)
                ax.text(self.spot["Row"].loc[index] - 1.25, self.spot["Column"].loc[index] - 1.25, self.spot["Virus"].loc[index],
                        va='center', ha='center', fontsize=8, weight="bold", rotation=45)

        # set tick marks for grid
        ax.set_xticks(self.spot["Row"].values - 1.5)
        ax.set_yticks(self.spot["Column"].values - 1.5)

        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_xlim(x.min() - 0.5, x.max() + 0.5)
        ax.set_ylim(y.max() + 0.5, y.min() - 0.5)
        ax.grid()

        return fig

    def barplot(self, align="pep", **kwargs):
        """ Create barplot

        :param align: flag for aligning with peptides 'pep' or virus 'vir'
        :param kwargs:
        :return:
        """
        fig = plt.figure(**kwargs)
        ax2 = plt.subplot(111)
        ax = ax2.twiny()

        # for x-axis ticks and labels
        peptide_ticks_x_axis = []
        peptide_label_x_axis = []

        virus_ticks_x_axis = []
        virus_label_x_axis = []

        cmap = plt.get_cmap('jet')
        Nvirus = len(self.spot["Virus"].unique())
        Npeptides = len(self.spot["Peptide"].unique())

        # peptide aligned
        if align == "pep":
            spacing = 1.0/Nvirus
            for index_peptide, peptide in enumerate(self.spot["Peptide"].unique()):
                for index_virus, virus in enumerate(self.spot["Virus"].unique()):
                    plt.axvline(x=[index_peptide + index_virus * spacing + spacing * 0.5 ], linestyle='dashed')
                    # add x-tick position and label for virus
                    virus_ticks_x_axis.append(index_peptide + index_virus * spacing)
                    virus_label_x_axis.append(virus)
                    data=self.spot.loc[(self.spot['Peptide'] == peptide) & (self.spot["Virus"] == virus)]["Intensity"]
                    data_std=self.spot.loc[(self.spot['Peptide'] == peptide) & (self.spot["Virus"] == virus)]["Std"]

                    ax.scatter(index_peptide  * np.ones(data.shape) + index_virus * spacing, data,
                               s=500/(Nvirus*Npeptides), color='k', marker='o')
                    if len(data) > 1:
                        bp = ax.boxplot(data.values, positions=[index_peptide+index_virus*spacing],
                                         patch_artist=True, showfliers=False,widths=spacing)
                        plt.setp(bp['boxes'], color=cmap(index_virus*1.0 / Npeptides), alpha=0.5)

                # add x-tick position and label for peptide
                peptide_ticks_x_axis.append(index_peptide + 0.5 - spacing * 0.5)
                peptide_label_x_axis.append(peptide)
                #draw vertical line to seperate peptides
                plt.axvline(x=index_peptide + 1 - spacing * 0.5 )

            # setup upper x-axis
            plt.xticks(peptide_ticks_x_axis,peptide_label_x_axis,fontsize="large",rotation=90)
            ax.set_xlim(-0.5*spacing , index_peptide + 1 - spacing * 0.5 )

            # setup lower x-axis
            ax2.set_xlim(ax.get_xlim())
            ax2.set_xticks(virus_ticks_x_axis)
            #["peptide " + str(s) for s in output.columns]
            ax2.set_xticklabels(virus_label_x_axis,rotation=90,fontsize="large")

        # virus aligned
        elif align == "vir":
            spacing = 1.0 / Npeptides
            for index_virus, virus in enumerate(self.spot["Virus"].unique()):
                for index_peptide, peptide in enumerate(self.spot["Peptide"].unique()):
                    plt.axvline(x=[index_virus + index_peptide * spacing + spacing * 0.5], linestyle='dashed')
                    # add x-tick position and label for virus
                    peptide_ticks_x_axis.append(index_virus + index_peptide * spacing)
                    peptide_label_x_axis.append(peptide)
                    data = self.spot.loc[(self.spot['Peptide'] == peptide) & (self.spot["Virus"] == virus)]["Intensity"]
                    ax.scatter(index_virus * np.ones(data.shape) + index_peptide * spacing, data,
                               s=500 / (Nvirus * Npeptides), color='k', marker='o')

                    if len(data) > 1:
                        bp = ax.boxplot(data.values, positions=[index_virus + index_peptide * spacing],
                                        patch_artist=True, showfliers=False, widths=spacing)
                        plt.setp(bp['boxes'], color=cmap(index_peptide * 1.0 / Npeptides), alpha=0.5)

                # add x-tick position and label for peptide
                virus_ticks_x_axis.append(index_virus + 0.5 - spacing * 0.5)
                virus_label_x_axis.append(virus)
                # draw vertical line to seperate peptides
                plt.axvline(x=index_virus + 1 - spacing * 0.5)

            # setup upper x-axis
            plt.xticks(virus_ticks_x_axis, virus_label_x_axis, fontsize="large", rotation=90)
            ax.set_xlim(-0.5 * spacing, index_virus + 1 - spacing * 0.5)

            # setup lower x-axis
            ax2.set_xlim(ax.get_xlim())
            ax2.set_xticks(peptide_ticks_x_axis)
            ax2.set_xticklabels(peptide_label_x_axis, rotation=90, fontsize="large")
            ax2.tick_params(labelsize="large")

        ax2.tick_params(labelsize="large")
        ax2.set_ylabel("Intensity", fontsize="x-large")

        return fig


    ### PCA ###

    def reshape_pca(self):
        # this should be later on in classfication
        spots_pca = self.spot[["Peptide", "Virus", "Intensity", "Replica","Referenz"]]
        # exclude referenz peptides
        spots_pca = spots_pca[~spots_pca["Referenz"]]
        spots_pca = spots_pca.set_index(["Virus", "Peptide", "Replica"])
        spots_pca = spots_pca.pivot_table(index=["Virus", "Replica"], columns="Peptide", values="Intensity")
        #  here a different method could be better. e.g. mean ...
        spots_pca_complete = spots_pca.dropna()
        # ---------------------------------------#
        # maybe an algorithm for selcting only
        # good peptides is reasonable

        return  spots_pca_complete

    def pca(self):
        """

        :return: (output_pca, pca_score, pca_components[0], feature_max_cor, cor_to_strongest_feature): where:
            output_pca: pd.Dataframe with output_pca[0,1] pca values, output_pca['color'] in numeric values
            pca_score: PCA score for principal components (explained variance)
            pca_components[0]: PCA Components in feature space
            feature_max_cor: Feature with max. Variation
            cor_to_strongest_feature: The Correlation of the Features to the Feature with highest Varaiance

        """
        # this should be probably later on in classification
        spots_pca_complete = self.reshape_pca()
        if len(spots_pca_complete.columns) > 1 and len(spots_pca_complete.index) > 1:
            # build model
            pca = decomposition.PCA()
            pca.fit(spots_pca_complete)

            # perform pca
            pca_score = pca.explained_variance_ratio_
            pca_components = pca.components_
            x_train_pca = pca.transform(spots_pca_complete)
            output_pca = pd.DataFrame(x_train_pca, index=spots_pca_complete.index)

            # store color information
            output_pca.insert(0, "color", utils.map_strings_to_number(spots_pca_complete.index.droplevel(1)))

            corr_values = dict(zip(spots_pca_complete.columns, pca_components[0]))
            feature_max_cor = max(corr_values, key=corr_values.get)
            cor_to_strongest_feature = spots_pca_complete[spots_pca_complete.columns[:]].corr()[feature_max_cor][:]


        elif len(spots_pca_complete.columns) < 2:
            raise Exception("You have selected measuraments with less than two peptides. No PCA is possible ")
        elif len(spots_pca_complete.index) < 2:
            raise Exception("You have selected measuraments with less than two two viruses (also replica) with complete set." \
                  "PCA possible  ")

        return output_pca, pca_score, pca_components[0], feature_max_cor, cor_to_strongest_feature

    ### CORRELATION ###

    def correlation_plot(self, type="both",**kwargs):
        spots_cor = self.spot[["Peptide", "Virus", "Intensity"]]
        spots_cor = spots_cor.set_index(["Virus", "Peptide"])
        unique_virus_peptide = spots_cor.pivot_table(index="Virus", columns="Peptide", values="Intensity")

        if type == "both":

            if unique_virus_peptide.shape[0] < 2 or unique_virus_peptide.shape[1] < 2:
                raise Exception("-> Not enough viruses or peptides in data set")
            else:
                data = unique_virus_peptide.T.corr()
                fig, (ax, ax2) = plt.subplots(1, 2, **kwargs)
                correlation.plot_corr_ellipses(data, ax=ax, cmap='seismic', clim=[-1, 1])
                ax.margins(0.1)
                m = correlation.plot_corr_ellipses(unique_virus_peptide.corr(), ax=ax2, cmap='seismic', clim=[-1, 1])
                cb = fig.colorbar(m)
                cb.ax.tick_params(labelsize="xx-large")
                cb.set_label('Correlation coefficient', fontsize="xx-large")
                ax2.margins(0.1)
                fig.tight_layout()

        elif type == "pep":
            if unique_virus_peptide.shape[1] < 2 :
                raise Exception("-> Not enough peptides in data set")
            else:
                fig, ax = plt.subplots( **kwargs)
                m = correlation.plot_corr_ellipses(unique_virus_peptide.corr(), ax=ax, cmap='seismic', clim=[-1, 1])
                cb = fig.colorbar(m)
                cb.ax.tick_params(labelsize="xx-large")
                cb.set_label('Correlation coefficient', fontsize="xx-large")
                ax.margins(0.1)
                fig.tight_layout()
        elif type == "vir":
            if unique_virus_peptide.shape[0] < 2:
                raise Exception("-> Not enough viruses in data set")
            else:
                fig, ax = plt.subplots(**kwargs)
                m = correlation.plot_corr_ellipses(unique_virus_peptide.T.corr(), ax=ax, cmap='seismic', clim=[-1, 1])
                cb = fig.colorbar(m)
                cb.ax.tick_params(labelsize="xx-large")
                cb.set_label('Correlation coefficient', fontsize="xx-large")
                ax.margins(0.1)
                fig.tight_layout()

        return fig

