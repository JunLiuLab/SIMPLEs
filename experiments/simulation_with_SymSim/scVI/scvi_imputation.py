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
def scvi_impute(seed: int = 1, platform: str = "umi") -> None:
    fnm: str = f"sim_{ncell}_{ngene}_{seed}_{platform}_.csv"
    save_path: PosixPath = here('./scVI/data/symsim')
    # fullpath:PosixPath = here('./scVI/data/symsim').joinpath(fnm)

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

    out_path = here("./simutool/jobs/scvi_result").joinpath(fnm)
    np.savetxt(out_path, impute_values, delimiter=",")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="scVI")
    parser.add_argument("--platform",
                        type=str,
                        default="umi",
                        help="umi or nonumi")
    parser.add_argument("--seed",
                        type=int,
                        default=1,
                        help="seed from 1 to 20")
    args = parser.parse_args()
    scvi_impute(args.seed, args.platform)
