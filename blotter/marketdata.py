import os
import pandas as pd


class MarketData():
    """
    Class for managing the data access layer. Provides functionality for
    reading data from disk and access methods
    """

    def __init__(self, **data_paths):
        """
        Parameters:
        -----------

        data_paths:
            key value pairs where the key is translated into the object
            attribute and the value is a string representing the path to
            the folder or file of data for all market data to be read.
            The attribute is either a dictionary with the file names as keys
            when a folder path is given or a pandas.DataFrame when a csv file
            path is given. The folder and all files in the folder should have
            the following structure

                folder/
                  instrument1.csv
                  instrument2.csv
                  ...

            where each csv file has the following structure

                datetime,CLZ6 Comdty
                2015-10-15T09:30:00, 35.76
                2015-10-15T09:35:00, 35.80
                2015-10-15T09:40:00, 35.70

            Note that each instrument csv file can contain one or more columns
            of market data

        Examples:
        ---------
        >>> mdata = MarketData(prices='data/prices', rates='data/rates')
        """
        for key in data_paths:
            path = data_paths[key]
            if os.path.isfile(path):
                data = self._read_file(path)
            else:
                files = os.listdir(path)
                fullfiles = [os.path.join(path, file) for file in files]
                data = dict()
                for file in fullfiles:
                    df = self._read_file(file)
                    name = os.path.split(file)[-1].split('.')[0]
                    data[name] = df
            setattr(self, key, data)

    @staticmethod
    def _read_file(file):
        return pd.read_csv(file, index_col=0, parse_dates=True)
