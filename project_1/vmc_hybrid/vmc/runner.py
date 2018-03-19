import tqdm
import pandas as pd
import numpy as np
from .resampling_methods import bootstrap, blocking

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

def run_all(sampler, parameters, parameter_names, bootstrap_samples, time=True,
        **sampler_kwargs):

    _quantities = ["energy", "variance", "std", "acceptance"]
    _boot_quantities = ["boot_var", "boot_std"]
    _block_quantities = ["block_var", "block_std"]

    columns = \
            parameter_names + _quantities + _boot_quantities + _block_quantities

    df = pd.DataFrame(columns=columns)

    for i, parameter in enumerate(parameter_names):
        df[parameter] = parameters[:, i]

    df = df.astype("double")

    for i in tqdm.tqdm(range(len(parameters))) \
            if time else range(len(parameters)):

        # Redistribute and set variational parameters in wavefunction
        sampler.initialize(df.iloc[i][parameter_names].values)

        sampler.sample(**sampler_kwargs)

        df.loc[i, _quantities] = [
            sampler.get_energy(), sampler.get_variance(), sampler.get_stddev(),
            sampler.get_acceptance_ratio()]

        local_energies = sampler.get_local_energies()

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
