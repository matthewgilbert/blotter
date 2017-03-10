import pandas as pd


def price_date_fill(dates, input_files, output_files):
    """
    Reindex a set of csv files to a new DatetimeIndex and fill forward values

    Parameters:
    -----------
    dates: pd.DatetimeIndex
        Dates to reindex prices to
    input_files: list
        List of strs of csv file names to reindex
    output_files: list
            List of strs of csv file names to write reindexed prices
    """
    for fi, fo in zip(input_files, output_files):
        print("Processessing file: %s" % fi)
        df = pd.read_csv(fi, parse_dates=True, index_col=0)
        idx_name = df.index.name
        df = df.reindex(dates).fillna(method="ffill")
        df.index.name = idx_name
        print("Writing file: %s" % fo)
        df.reset_index().to_csv(fo, index=False)
