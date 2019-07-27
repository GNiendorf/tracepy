import pandas as pd

def save_optics(geo_params, filename):
    """ Save geometry to an optics table in csv format. """
    table = pd.DataFrame.from_dict(geo_params)
    table[['X', 'Y', 'Z']] = pd.DataFrame(table.P.values.tolist(), index= table.index)
    table[['Alpha', 'Beta', 'Gamma']] = pd.DataFrame(table.D.values.tolist(), index= table.index)
    del[table['P']]; del[table['D']]
    pd.DataFrame.to_csv(table, filename)
