import pandas as pd

# TODO: Rewrite save_optics to iterate over geo_params to handle edge
#       cases, and convert to dataframe at end before exporting to csv.
def save_optics(geo_params, filename: str):
    """Save geometry to an optics table in csv format.

    Note
    ----
    Directly saves the csv file rather than returning something.

    Parameters
    ----------
    geo_params : list of dictionaries
        Dictionaries correspond to surfaces.
    filename : str
        Output filename.

    """

    table = pd.DataFrame.from_dict(geo_params)
    table[['X', 'Y', 'Z']] = pd.DataFrame(table.P.values.tolist(), index= table.index)
    table[['Alpha', 'Beta', 'Gamma']] = pd.DataFrame(table.D.values.tolist(), index= table.index)
    del[table['P']]; del[table['D']]
    pd.DataFrame.to_csv(table, filename)
