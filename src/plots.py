import matplotlib.pyplot as plt

class Plots:

    def __init__(self):
        pass

    def make_genome_ideogram(self, AF, genome_baf, Dosage, capture):
        '''  
            Plot dosage,baf ideograms for the whole genome and return the plt object
            Allows you to dump in a PDF file or save as an image 
        '''
        n_rows, n_cols = 2, 1
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(24, 12), layout="compressed", sharex=True)
        axes[0].set_ylim([-0.1, 2.1])
        axes[0].set_yticks(np.arange(0, 2, 0.25))
        Dosage.plot_ideogram_ax(axes[0], capture)
        # Set custom Y-axis ticks
        custom_ticks = [0, 0.25, 0.337, 0.5, 0.667, 0.75, 1]
        axes[1].set_yticks(custom_ticks, labels=[str(tick) for tick in custom_ticks])
        AF.plot_baf_ideogram(axes[1], genome_baf, capture)
        plt.tight_layout()

        return plt