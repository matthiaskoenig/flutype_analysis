import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.collections import EllipseCollection

from sklearn import decomposition

import os


def load_data(data_id):

    """
    :param data_id: The ID of the Data
    :return: gal_pep,gal_vir,intensity,meta
    """

    print("-" * 80)
    print("Loading data corresponding to data ID :{}".format(data_id))
    print("-" * 80)

    count_vir = 0
    count_pep = 0
    count_meta = 0
    count_int = 0
    count_std = 0

    for fname in os.listdir('data/'):  # change directory as needed
        if os.path.isfile('data/{}'.format(fname)):  # make sure it's a file,  not a directory entry
            if "{}_vir.gal".format(data_id) in fname:  # search for vir.gal
                gal_vir = pd.read_csv("data/{}".format(fname),
                                      sep='\t', index_col="ID")

                print("Virus .gal :{}".format(fname))
                count_vir += 1

            elif "{}_pep.gal".format(data_id) in fname: # search for pep.gal
                gal_pep = pd.read_csv("data/{}".format(fname), sep='\t',
                                      index_col="ID")
                print("Peptide .gal :{}".format(fname))
                count_pep += 1

            elif "{}.meta".format(data_id) in fname: # search for meta
                meta = pd.read_csv("data/{}".format(fname), sep='\t')

                print("Meta  :{}".format(fname))
                count_meta += 1

            elif "{}.csv".format(data_id) in fname: # search for intensity
                spot_intensity = pd.read_csv("data/{}".format(fname),
                                             sep='\t', index_col=0)
                print("Spot intensity file  :{}".format(fname))
                count_int += 1

            elif "{}_std.csv".format(data_id) in fname:  # search for intensity
                spot_intensity_std = pd.read_csv("data/{}".format(fname),
                                             sep='\t', index_col=0)
                print("Intensity standard deviation:{}".format(fname))
                count_std += 1

    counter_all = [count_pep,count_meta,count_vir,count_int]

    print("-" * 80)
    if not (count_pep == 1):
        if count_pep < 1:
            raise Exception(
                "The pep.gal data was not found in /data")
        elif count_pep > 1:
            raise Exception(
            "To many coresponding pep.gal data in /data")


    if not (count_vir == 1):

        if count_vir < 1:
            raise Exception(
                "The vir.gal data was not found in /data")
        elif count_vir > 1:
            raise Exception(
                "To many coresponding vir.gal data in /data.")

    if not (count_meta == 1):

        if count_meta < 1:
            raise Exception(
                "The meta data was not found in /data")
        elif count_meta > 1:
            raise Exception(
                "To many coresponding meta data in /data.")

    if not (count_int == 1):
        if count_int < 1:
            raise Exception(
                "The intensity data was not found in /data")

        elif count_int > 1:
            raise Exception(
                "To many coresponding data with Spot intensity in /data.")
    elif not (count_std == 1):
            if count_std < 1:
                spot_intensity_std = 0
                print("Spot intensities for the data ID ({}) are not averaged but primary values".format(data_id))
                print("-" * 80)
            elif count_std > 1:
                raise Exception(
                    "To many coresponding data with standard deviation in /data.")

    if all( count == 1 for count in  counter_all):
        print(" Necessary coresponding data was loaded ")
        print("-" * 80)
        return [data_id,spot_intensity,spot_intensity_std ,gal_pep ,gal_vir ,meta]


def map_strings_to_number(strings):
    """Transforms list of strings into numbers."""
    counter = 1
    d = {}
    numbers = []
    for s in strings:
        if s in d:
            pass
        else:
            d[s] = counter
            counter+=1
        numbers.append(d[s])
    return numbers


def discrete_cmap(N, base_cmap=None):
    """Create an N-bin discrete colormap from the specified input map"""
    # Note that if base_cmap is a String or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for String, None, or a colormap instance:
    if N == 1:
        base = plt.cm.get_cmap(base_cmap)
        color_list = base(np.linspace(0, 1, N + 1))
        cmap_name = base.name + str(N)

    else:

        base = plt.cm.get_cmap(base_cmap)
        color_list = base(np.linspace(0, 1, N))
        cmap_name = base.name + str(N)

    return base.from_list(cmap_name, color_list, N)

def pca_2dplot(output_pca,**kwargs):
        fig = plt.figure(**kwargs)
        plt.xlabel("princple component 1", size='large')
        plt.ylabel("princple component 2", size='large')
        plt.tick_params(labelsize="large")
        plt.scatter(output_pca[0], output_pca[1],
                    c=output_pca["color"].values, cmap=discrete_cmap(output_pca["color"].values.max(), 'jet'))

        color_mapping = output_pca[['color']].sort_values(by=["color"])

        cbar = plt.colorbar(ticks=sorted(output_pca["color"].unique()))
        cbar.set_ticklabels(color_mapping.index.droplevel(1).unique())
        cbar.ax.tick_params(labelsize="large")
        plt.clim(0.5, output_pca["color"].values.max() + 0.5)

        return fig

def plot_corr_ellipses(data, ax=None, **kwargs):

    M = np.array(data)
    if not M.ndim == 2:
        raise ValueError('data must be a 2D array')
    if ax is None:
        fig, ax = plt.subplots(1, 1, subplot_kw={'aspect':'equal'})
        ax.set_xlim(-0.5, M.shape[1] - 0.5)
        ax.set_ylim(-0.5, M.shape[0] - 0.5)

    # xy locations of each ellipse center
    xy = np.indices(M.shape)[::-1].reshape(2, -1).T

    # set the relative sizes of the major/minor axes according to the strength of
    # the positive/negative correlation
    w = np.ones_like(M).ravel()
    h = 1 - np.abs(M).ravel()
    a = 45 * np.sign(M).ravel()

    ec = EllipseCollection(widths=w, heights=h, angles=a, units='x', offsets=xy,
                           transOffset=ax.transData, array=M.ravel(), **kwargs)
    ax.add_collection(ec)

    # if data is a DataFrame, use the row/column names as tick labels
    if isinstance(data, pd.DataFrame):
        ax.set_xticks(np.arange(M.shape[1]))
        ax.set_xticklabels(data.columns, rotation=90,fontsize="xx-large")
        ax.set_yticks(np.arange(M.shape[0]))
        ax.set_yticklabels(data.index,fontsize="xx-large")

    return ec

class Analysis():

    def __init__(self, data_id, spot_intensity,spot_intensity_std ,gal_pep ,gal_vir ,meta):

        '''
        :param gal_pep:
        :param gal_vir:
        :param intensity:
        :param meta:
        '''

        self.gal_pep = gal_pep
        self.gal_vir = gal_vir
        self.spot_intensity = spot_intensity
        self.spot_intensity_std = spot_intensity_std
        self.meta = meta
        self.data_id = data_id
        self.spot = self.spot()

    def spot(self):
        vir_cor = self.gal_vir.pivot(index="Row", columns="Column", values="Name")
        pep_cor = self.gal_pep.pivot(index="Row", columns="Column", values="Name")

        # merge complete spotinformation
        vir_cor_unstacked = vir_cor.unstack()
        spot = pep_cor.unstack()
        spot = spot.reset_index()
        spot = spot.rename(columns={0: "Peptide"})
        spot["Referenz"] = spot["Peptide"].str.contains("Leuchte")
        spot["Virus"] = vir_cor_unstacked.values
        spot["Intensity"] = self.spot_intensity.unstack().values

        if type(self.spot_intensity_std) is int:
            spot["Std"] = self.spot_intensity_std
        else:
            spot["Std"] = self.spot_intensity_std.unstack().values

        spot["Replica"] = 0

        for virus_unique in spot["Virus"].unique():
            for peptide_unique in spot["Peptide"].unique():
                replica = 0
                for index in spot.index:
                    if spot["Virus"][index] == virus_unique and spot["Peptide"][index] == peptide_unique:
                        spot.set_value(index, "Replica", replica)
                        replica += 1




        return spot

    def heatmap(self, descript = True , heatmap = True , **kwargs):

        fig, ax = plt.subplots(**kwargs)
        # imshow portion
        if descript:
            alpha=0.5
        else:
            alpha=1.0
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

    def barplot(self, align="pep" ,**kwargs):

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
            output_pca.insert(0, "color", map_strings_to_number(spots_pca_complete.index.droplevel(1)))

            corr_values = dict(zip(spots_pca_complete.columns, pca_components[0]))
            feature_max_cor = max(corr_values, key=corr_values.get)
            cor_to_strongest_feature = spots_pca_complete[spots_pca_complete.columns[:]].corr()[feature_max_cor][:]


        elif len(spots_pca_complete.columns) < 2:
            raise Exception("You have selected measuraments with less than two peptides. No PCA is possible ")
        elif len(spots_pca_complete.index) < 2:
            raise Exception("You have selected measuraments with less than two two viruses (also replica) with complete set." \
                  "PCA possible  ")

        return output_pca, pca_score, pca_components[0], feature_max_cor, cor_to_strongest_feature


    def correlation_plot(self,type="both",**kwargs):
        spots_cor = self.spot[["Peptide", "Virus", "Intensity"]]
        spots_cor = spots_cor.set_index(["Virus", "Peptide"])
        unique_virus_peptide = spots_cor.pivot_table(index="Virus", columns="Peptide", values="Intensity")

        if type == "both":

            if unique_virus_peptide.shape[0] < 2 or unique_virus_peptide.shape[1] < 2:
                raise Exception("-> Not enough viruses or peptides in data set")


            else:
                data = unique_virus_peptide.T.corr()
                fig, (ax, ax2) = plt.subplots(1, 2, **kwargs)
                plot_corr_ellipses(data, ax=ax, cmap='seismic', clim=[-1, 1])
                ax.margins(0.1)
                m = plot_corr_ellipses(unique_virus_peptide.corr(), ax=ax2, cmap='seismic', clim=[-1, 1])
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
                m = plot_corr_ellipses(unique_virus_peptide.corr(), ax=ax, cmap='seismic', clim=[-1, 1])
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
                m = plot_corr_ellipses(unique_virus_peptide.T.corr(), ax=ax, cmap='seismic', clim=[-1, 1])
                cb = fig.colorbar(m)
                cb.ax.tick_params(labelsize="xx-large")
                cb.set_label('Correlation coefficient', fontsize="xx-large")
                ax.margins(0.1)
                fig.tight_layout()


        return fig

