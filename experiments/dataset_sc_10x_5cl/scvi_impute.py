import argparse
import torch
import numpy as np
import pandas as pd
from pyprojroot import here
from pathlib import PosixPath

from scvi.dataset import CsvDataset
from scvi.models import VAE
from scvi import set_seed
from scvi.inference import UnsupervisedTrainer, load_posterior

# * commons

set_seed(0)

# ** data
ncell = 300
ngene = 1000

# ** model
n_epochs = 400
lr = 1e-3
use_cuda = True


# * impute
def scvi_impute() -> None:
    fnm: str = "sc_10x_5cl_forimput_cnt.csv"
    save_path: PosixPath = here('./10xGenomics/scRNAseq')

    symsim_dataset = CsvDataset(fnm, save_path=save_path, gene_by_cell=True)

    vae = VAE(symsim_dataset.nb_genes)

    trainer = UnsupervisedTrainer(vae,
                                  symsim_dataset,
                                  train_size=1.0,
                                  use_cuda=use_cuda,
                                  frequency=5)

    trainer.train(n_epochs=n_epochs, lr=lr)

    full = trainer.create_posterior(trainer.model,
                                    symsim_dataset,
                                    indices=np.arange(len(symsim_dataset)))
    impute_values = full.sequential().imputation()

    outfnm: str = "scvi_impt.csv"
    out_path = here("./10xGenomics/impt/").joinpath(outfnm)
    np.savetxt(out_path, impute_values, delimiter=",")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="scVI")
    scvi_impute()
