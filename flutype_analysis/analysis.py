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
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from sklearn import decomposition
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

from flutype_analysis import correlation
from flutype_analysis import utils
from flutype_analysis import base

from matplotlib.figure import Figure

# from jinja2 import Environment, FileSystemLoader

import os


class Analysis(base.Base):

    def __init__(self, d_data):
        base.Base.__init__(self, d_data)


    def heatmap(self, descript=True, heatmap=True, **kwargs):
        """ Creates heatmap

        :param descript: boolean flag to show peptide information
        :param heatmap: boolean flag to show heatmap
        :param kwargs:
        :return:
        """
        fig = Figure(**kwargs)
        ax = fig.add_subplot(111)

        # imshow portion
        if descript:
            alpha = 0.5
        else:
            alpha = 1.0

        if heatmap:
            ax.imshow(self.spot.pivot(index='Row', columns='Column',
                                      values='Intensity'), interpolation='nearest', cmap="hot", alpha=alpha)
        # plt.pcolor(Spot["Intensity"].unstack())
        # text portion
        x, y = np.meshgrid(self.spot["Column"].unique() - 1, self.spot["Row"].unique() - 1)
        for index in self.spot.index:
            # print(y_val)
            if descript:
                ax.text(self.spot["Column"].loc[index] - 0.90, self.spot["Row"].loc[index] - 0.90, self.spot["Peptide"].loc[index],
                        va='center', ha='center', fontsize=8, weight="bold", rotation=45)
                ax.text(self.spot["Column"].loc[index] - 1.25, self.spot["Row"].loc[index] - 1.25, self.spot["Virus"].loc[index],
                        va='center', ha='center', fontsize=8, weight="bold", rotation=45)

        # set tick marks for grid
        ax.set_xticks(self.spot["Column"].values - 1.5)
        ax.set_yticks(self.spot["Row"].values - 1.5)

        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_xlim(x.min() - 0.5, x.max() + 0.5)
        ax.set_ylim(y.max() + 0.5, y.min() - 0.5)
        ax.grid()

        return fig

    def barplot(self, align="pep", scale="log", **kwargs):
        """ Create barplot

        :param align: flag for aligning with peptides 'pep' or virus 'vir'
        :param kwargs:
        :return:
        """
        fig = Figure(**kwargs)
        ax2 = fig.add_subplot(111)
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
                        plt.setp(bp['boxes'], color=cmap(index_virus*1.0 / Nvirus), alpha=0.5)

                # add x-tick position and label for peptide
                peptide_ticks_x_axis.append(index_peptide + 0.5 - spacing * 0.5)
                peptide_label_x_axis.append(peptide)
                #draw vertical line to seperate peptides
                plt.axvline(x=index_peptide + 1 - spacing * 0.5 )

            # setup upper x-axis
            plt.xticks(peptide_ticks_x_axis,peptide_label_x_axis,fontsize="x-large",rotation=90)
            ax.set_xlim(-0.5*spacing , index_peptide + 1 - spacing * 0.5 )


            # setup lower x-axis
            ax2.set_xticks(virus_ticks_x_axis)
            #["peptide " + str(s) for s in output.columns]
            ax2.set_xticklabels(virus_label_x_axis,rotation=90, fontsize="x-small" )

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
            plt.xticks(virus_ticks_x_axis, virus_label_x_axis, fontsize="x-large", rotation=90)
            ax.set_xlim(-0.5 * spacing, index_virus + 1 - spacing * 0.5)

            # setup lower x-axis
            ax2.set_xticks(peptide_ticks_x_axis)
            ax2.set_xticklabels(peptide_label_x_axis, rotation=90, fontsize = "x-small")

        ax.set_ylim(self.spot["Intensity"].min(), self.spot["Intensity"].max())
        ax2.set_yscale(scale)
        ax2.set_xlim(ax.get_xlim())
        ax2.set_ylabel("Intensity", fontsize="x-large")

        return fig


    ### PCA ###

    def reshape_pca(self, antibody=False):
        # this should be later on in classfication
        spots_pca = self.spot[["Peptide", "Virus", "Intensity", "Replica","Reference","No_Peptide","Buffer"]]
        # exclude referenz peptides
        spots_pca = spots_pca[~spots_pca["Reference"]]
        #exclude empty
        spots_pca = spots_pca[~spots_pca["No_Peptide"]]
        #exclude buffer
        spots_pca = spots_pca[~spots_pca["Buffer"]]

        #exclude antibodies
        if not antibody:
            spots_pca = spots_pca[~spots_pca["Peptide"].str.contains("AK")]


        spots_pca = spots_pca.set_index(["Virus", "Peptide", "Replica"])
        spots_pca = spots_pca.pivot_table(index=["Virus", "Replica"], columns="Peptide", values="Intensity")
        #  here a different method could be better. e.g. mean ...
        spots_pca_complete = spots_pca.dropna()


        # ---------------------------------------#
        # maybe an algorithm for selcting only
        # good peptides is reasonable



        return  spots_pca_complete

    def classifier(self, type ="PCA", **kwargs):
        """

        type    PCA :
                IPA : Independened Component Analysis
                LDA : LinearDiscriminantAnalysis

        :return: (output_pca, pca_score, pca_components[0], feature_max_cor, cor_to_strongest_feature): where:
            output_pca: pd.Dataframe with output_pca[0,1] pca values, output_pca['color'] in numeric values
            pca_score: PCA score for principal components (explained variance)
            pca_components[0]: PCA Components in feature space
            feature_max_cor: Feature with max. Variation
            cor_to_strongest_feature: The Correlation of the Features to the Feature with highest Varaiance

        """
        # this should be probably later on in classification
        spots_pca_complete = self.reshape_pca(**kwargs)

        if len(spots_pca_complete.columns) > 1 and len(spots_pca_complete.index) > 1:
            # build model
            if type == "PCA":
                classifier = decomposition.PCA()
                classifier.fit(spots_pca_complete)
            elif type == "ICA":
                classifier = decomposition.FastICA()
                classifier.fit(spots_pca_complete)

            elif type == "LDA":
                classifier = LinearDiscriminantAnalysis(n_components=5)
                y=spots_pca_complete.index.droplevel(1).values
                classifier.fit(spots_pca_complete,y)




            x_train_pca = classifier.transform(spots_pca_complete)
            output_pca = pd.DataFrame(x_train_pca, index=spots_pca_complete.index)

            # store color information
            output_pca.insert(0, "color", utils.map_strings_to_number(spots_pca_complete.index.droplevel(1)))




            if type == "ICA":
                pca_score = None


            else:
                pca_score = classifier.explained_variance_ratio_


            if type == "LDA":
                pca_components = 0
                pca_components = classifier.coef_
                two_pca_components = pd.DataFrame(pca_components[:2, :].T, index=spots_pca_complete.columns.values,
                                                  columns=["Pca_Component_1", "Pca_Component_2"])
                two_pca_components.sort_values(["Pca_Component_1"], inplace=True, ascending=False)

                return output_pca, pca_score,two_pca_components

            else:
                pca_components = classifier.components_
                two_pca_components = pd.DataFrame(pca_components[:2,:].T,index=spots_pca_complete.columns.values, columns=["Pca_Component_1","Pca_Component_2"])
                two_pca_components.sort_values(["Pca_Component_1"], inplace=True, ascending=False)
                corr_values = dict(zip(spots_pca_complete.columns, pca_components[0]))
                feature_max_cor = max(corr_values, key=corr_values.get)
                cor_to_strongest_feature = spots_pca_complete[spots_pca_complete.columns[:]].corr()[feature_max_cor][:]

                return output_pca, pca_score, two_pca_components, feature_max_cor, cor_to_strongest_feature




        elif len(spots_pca_complete.columns) < 2:
            raise Exception("You have selected measuraments with less than two peptides. No PCA is possible ")
        elif len(spots_pca_complete.index) < 2:
            raise Exception("You have selected measuraments with less than two two viruses (also replica) with complete set." \
                  "PCA possible  ")



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



    def components_to_html(self, **kwargs):
        """
        writes a pdf with with the principle components into the corresponding results directory.
        :return:
        """

        data_pca = self.classifier(**kwargs)
        env = Environment(loader=FileSystemLoader('..'))
        directory = "/html/"
        hfile = os.path.join(directory, "principle_component.html")
        template = env.get_template(hfile)
        template_vars = {"pca": data_pca[2].to_html()}
        html_out = template.render(template_vars)

        return html_out

