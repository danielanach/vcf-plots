import pandas as pd

class SubstitutionSpectrum:
    def __init__(self,samples):
        self.samples = samples
        self.mut_count_dict = {}

        muts = ['C→A',
                'C→G',
                'C→T',
                'T→A',
                'T→C',
                'T→G',]
        for mut in muts:
            self.mut_count_dict[mut] = []

    def normalized_muts(self):
        """Normalize substitution type counts to sum to 1 per sample.
        Parameters
        ----------
        vcf : str
            VCF file path.

        sample_list : list(str)
            list of sample names.

        n_bins : int
            number of bins to split variants by VAF into
        """

        mut_df = pd.DataFrame.from_dict(self.mut_count_dict)
        mut_df_norm = mut_df.divide(mut_df.sum(axis=1),axis='rows')
        mut_norm_dict = mut_df_norm.to_dict('list')

        return mut_norm_dict

    def counts_by_category(self):
            """Normalize substitution type counts to sum to 1 per sample.
            Parameters
            ----------

            """

            mut_df = pd.DataFrame.from_dict(self.mut_count_dict)
            return mut_df.sum(axis=1)
