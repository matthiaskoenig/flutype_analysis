import unittest


class MyTestCase(unittest.TestCase):

    """
    Test all functions in flutype-analysis
    """

    # ----------------------------------------------------------------- #
    # ----------------------quality_control.py------------------------- #
    # ----------------------------------------------------------------- #
    # ----------------shell inherit from image2numeric---------------- #

    def test_is_path_variable_occupied(self):
        """
        test for name conflicts
        :return:
        """
        self.assertEqual(True, False)

    def test_data_images_complete(self):
        self.assertEqual(True, False)

    def test_is_standard_format_image(self):
        self.assertEqual(True, False)

    def test_align_resolution(self):
        self.assertEqual(True, False)

    def test_align_format(self):
        self.assertEqual(True, False)

    def test_transform2standard_image(self):
        self.assertEqual(True, False)

    def test_are_all_fluor_spots_found(self):
        self.assertEqual(True, False)

    def test_are_all_spots_found(self):
        self.assertEqual(True, False)

    def test_quality_of_spot(self):
        self.assertEqual(True, False)

    def test_quality_of_grid(self):
        self.assertEqual(True, False)

    def test_std_of_duplicates(self):
        self.assertEqual(True, False)

    def test_quality_of_array(self):
        # result of quality of spot, grid and std of duplicates
        # add to meta
        self.assertEqual(True, False)

    # ------------------shell inherit from analysis------------------ #

    def test_is_standard_format_numeric(self):
        self.assertEqual(True, False)

    def test_transform2standard_numeric(self):
        self.assertEqual(True, False)

    def test_is_complementary_data_present(self):
        # meta data, pep.gal ,vir.gal
        self.assertEqual(True, False)

    def test_std_of_duplicates_mtp(self):
        self.assertEqual(True, False)

    def test_quality_of_mtp(self):
        # result of std of duplicates
        self.assertEqual(True, False)

    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # ----------------------image2numeric.py-------------------------- #
    # ---------------------------------------------------------------- #

    # load images before, after incubation and after washing:
    # Total 3 files per measurement.
    def test_load_data_images(self):
        self.assertEqual(True, False)

    def test_find_grid(self):
        #
        self.assertEqual(True, False)

    # for spot recognition for "Leuchtefix"
    def test_find_fluor_spots(self):
        self.assertEqual(True, False)

    # for spot recognition for all spots after incubation
    def test_find_all_spots(self):
        self.assertEqual(True, False)

    # after spot recognition
    def test_generate_data(self):
        self.assertEqual(True, False)

    # ---------------------------------------------------------------- #
    # -------------------------analysis.py---------------------------- #
    # ---------------------------------------------------------------- #

    # after data generation with spot recognition
    def test_load_data(self):
        self.assertEqual(True, False)

    def test_complete_spot_data(self):
    # box,column,row,intensity,intensity_std,pep,vir,vir_color,replica,Leuchtefix
        self.assertEqual(True, False)

    def test_heat_map(self):
        self.assertEqual(True, False)

    def test_heat_map_mini(self):
        self.assertEqual(True, False)

    def test_boxplot_vir(self):
        self.assertEqual(True, False)

    def test_boxplot_pep(self):
        self.assertEqual(True, False)

    def test_map_string2number(self):
        self.assertEqual(True, False)

    def test_discrete_cmap(self):
        self.assertEqual(True, False)

    def test_correlation_matrix_peptide(self):
        self.assertEqual(True, False)

    def test_correlation_matrix_virus(self):
        self.assertEqual(True, False)

    # ---------------------------------------------------------------- #
    # ----------------------classification.py----------- ------------- #
    # ---------------------------------------------------------------- #

    def test_vir_pep_model(self):
        self.assertEqual(True, False)

    def test_predict_virus(self):
        # virus, certainty
        self.assertEqual(True, False)

if __name__ == '__main__':
    unittest.main()
