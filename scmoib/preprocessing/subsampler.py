from sklearn.model_selection import train_test_split


class Subsampler:
    def __init__(self, adata, seed=0):
        self.df = adata.obs.copy()
        self.tmp_df = adata.obs.copy()
        self.history = []
        self.seed = seed

    def clear(self):
        self.df = self.tmp_df.copy()
        self.history = []
        print("Returned to the initial state")

    def show_history(self):
        print("\n".join(self.history))

    def split(self, omic_batch, cell_type, drop_size=0.2):
        omic_df = self.df[self.df['omic_batch'] == omic_batch]
        res = train_test_split(omic_df.index,
                               omic_df[cell_type],
                               test_size=drop_size,
                               random_state=self.seed,
                               shuffle=True,
                               stratify=omic_df[cell_type])[0]
        inds = list(res) + list(self.df[self.df['omic_batch'] != omic_batch].index)
        self.df = self.df.loc[inds, :]
        self.history.append(f"Remove {drop_size * 100}% of {omic_batch} cells")

    def filter_cell_type(self, omic_batch, cell_type, drop_type, drop_level=1):
        inds = list(self.df[(self.df['omic_batch'] == omic_batch) & (self.df[cell_type] != drop_type)].index)
        inds += list(self.df[self.df['omic_batch'] != omic_batch].index)
        if drop_level < 1:
            cell_df = self.df[(self.df['omic_batch'] == omic_batch) & (self.df[cell_type] == drop_type)]
            res = train_test_split(cell_df.index,
                                   cell_df[cell_type],
                                   test_size=drop_level,
                                   random_state=self.seed,
                                   shuffle=True,
                                   stratify=cell_df[cell_type])[0]
            inds += list(res)
        self.df = self.df.loc[inds, :]
        self.history.append(f"Remove {drop_level * 100}% of {drop_type} cells from {omic_batch}")

    def subsample(self, adata):
        inds = list(self.df.sort_values(by=['omic_batch', 'index']).index)
        return adata[inds].copy()
