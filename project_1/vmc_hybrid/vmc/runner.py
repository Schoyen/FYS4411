import tqdm
import pandas as pd

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
