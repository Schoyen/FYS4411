import tqdm
import time
import pandas as pd
import numpy as np
from .resampling_methods import bootstrap, blocking
import numba

def sample_local_energies(sampler, wavefunction, parameters, time=True,
        **sampler_kwargs):

    if not "sample_local_energies" in sampler_kwargs:
        sampler_kwargs["sample_local_energies"] = True

    # Use index into parameters array as column names
    df = pd.DataFrame(columns=[str(i) for i in range(parameters.size)])

    for i in tqdm.tqdm(range(parameters.size)) \
            if time else range(parameters.size):

        wavefunction.redistribute()
        wavefunction.set_parameters(parameters[i])

        sampler.sample(**sampler_kwargs)

        df[str(i)] = sampler.get_local_energies()

    return df

@numba.jit
def run_all(sampler, parameters, parameter_names, bootstrap_samples,
        timeit=True, **sampler_kwargs):

    if not "sample_local_energies" in sampler_kwargs:
        sampler_kwargs["sample_local_energies"] = True

    _quantities = ["energy", "variance", "std", "acceptance", "sampling_time"]

    _boot_quantities = ["boot_var", "boot_std"] if bootstrap_samples > 0 else []

    _block_quantities = ["block_var", "block_std"]

    columns = \
            parameter_names + _quantities + _boot_quantities + _block_quantities

    df = pd.DataFrame(columns=columns)

    for i, parameter in enumerate(parameter_names):
        df[parameter] = parameters[:, i]

    df = df.astype("double")

    for i in tqdm.tqdm(range(len(parameters))) \
            if timeit else range(len(parameters)):

        # Redistribute and set variational parameters in wavefunction
        sampler.initialize(df.iloc[i][parameter_names].values)

        t0 = time.time()
        sampler.sample(**sampler_kwargs)
        t1 = time.time()

        df.loc[i, _quantities] = [
            sampler.get_energy(), sampler.get_variance(), sampler.get_stddev(),
            sampler.get_acceptance_ratio(), t1 - t0]

        local_energies = sampler.get_local_energies()

        if bootstrap_samples > 0:
            df.loc[i, _boot_quantities] = bootstrap(
                    local_energies, bootstrap_samples)

        df.loc[i, _block_quantities] = blocking(local_energies)

    return df

def run_experiment(sampler, wavefunction, parameters, parameter_names,
        time=True, **sampler_kwargs):

    _quantities = ["energy", "variance", "acceptance"]

    df = pd.DataFrame(columns=parameter_names)

    for i, parameter in enumerate(parameter_names):
        df[parameter] = parameters[:, i]

    for quantity in _quantities:
        df[quantity] = 0

    for i in tqdm.tqdm(range(len(parameters))) \
            if time else range(len(parameters)):

        wavefunction.redistribute()
        wavefunction.set_parameters(df.iloc[i][parameter_names].values)

        sampler.sample(**sampler_kwargs)

        df.loc[i, _quantities] = [
            sampler.get_energy(), sampler.get_variance(),
            sampler.get_acceptance_ratio()]

    return df
