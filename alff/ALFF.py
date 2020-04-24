# -*- coding: utf-8 -*-
import os
import sys
import asyncio
import argparse
import requests
import nest_asyncio
import pandas as pd
from timeit import default_timer
from concurrent.futures import ThreadPoolExecutor

class ALFF:
    """ 
    Main class of Allele Frequency Finder.
    
    Contains variables from the command line and functions which use them.
    Running the programme with -h or --help will give a detailed breakdown of 
    arguments and their defaults.
    """

    def __init__(self):
        parser = argparse.ArgumentParser(description = "This script will call the NCBI variation API to append allele frequencies garnered from their ALFA project.\nFor more information, please see https://api.ncbi.nlm.nih.gov/variation/v0.")
        
        parser.add_argument("-i", 
                            "--input", 
                            type = str,
                            help = "Full or relative path to input file.", 
                            required = True, 
                            default = "")
        parser.add_argument("-is",
                            "--isep",
                            type = str,
                            help = "File separator for input file.",
                            required = False,
                            default = "\t")
        parser.add_argument("-o", 
                            "--output", 
                            help = "Full or relative path to output file. If the file exists, it will be overwritten.",
                            required = False,
                            default = "alff_output.txt")
        parser.add_argument("-os",
                            "--osep",
                            type = str,
                            help = "File separator for output file.",
                            required = False,
                            default = "\t")
        parser.add_argument("-snp",
                            "--snp_col",
                            type = str,
                            help = "Column header for SNP information. Will default to the first column if no string is given.",
                            required = False,
                            default = "")
        parser.add_argument("-allele",
                            "--allele_col",
                            type = str,
                            help = "Column header for allele information to find the frequency for. Will default to the second column if no string is given.",
                            required = False,
                            default = "")
        parser.add_argument("-org",
                            "--organism",
                            type = str,
                            help = "Organism for whose allele frequencies to search. Please see https://www.ncbi.nlm.nih.gov/bioproject/browse for more details. Defaults to human.",
                            required = False,
                            default = "PRJNA507278")
        parser.add_argument("-pop",
                            "--population",
                            type = str,
                            help = "Population subtype to use, defaults to European. Please see https://www.ncbi.nlm.nih.gov/snp/docs/gsr/alfa/ for more details.",
                            required = False,
                            default = "SAMN10492695")
        parser.add_argument("-w",
                            "--workers",
                            type = int,
                            help = "How many requests to send at once. Be mindful to not overload servers by setting too high. Setting this lower will take longer but is less likely to timeout.",
                            required = False,
                            default = 5)
        parser.add_argument("-t",
                            "--timeout",
                            type = int,
                            help = "(Seconds) How long to wait before a request is considered as timed out.",
                            required = False,
                            default = 5)
        parser.add_argument("-a",
                            "--attempts",
                            type = int,
                            help = "How many attempts should be tried if the request is met with a timeout or other related error.",
                            required = False,
                            default = 5)
        parser.add_argument("-v",
                            "--verbose",
                            type = int,
                            help = "(Boolean) Will print verbose warnings and messages. This will be automatically written to 'eaf.log' in the folder where this is run.",
                            required = False,
                            default = 0)
        
        argument = parser.parse_args()
        
        # New log file
        f = open("eaf.log", "w+")
        f.close()
        
        if argument.input:
            self.input = argument.input
        else:
            self.log("Error: Input file is required.", display=True)
            parser.print_help()
            sys.exit(1)
            
        self.isep = argument.isep
        self.osep = argument.osep
        self.output = argument.output
        self.snp_col = argument.snp_col
        self.allele_col = argument.allele_col
        
        self.organism = argument.organism
        self.population = argument.population
        
        self.workers = argument.workers
        self.timeout = argument.timeout
        self.attempts = argument.attempts
        self.verbose = argument.verbose
        
        self.d = {} # Temp holder for results
        self.run()
        
    def run(self):
        """
        Run the main body of Freq Finder.
        Will check if the SNP and Allele columns are usable and proceed to
        run the asynchronous API request loop.

        Returns
        -------
        None.

        """
        try:
            df = pd.read_csv(self.input, sep = self.isep)
        except FileNotFoundError:
            self.log("Error: File not found \"{}\". Please double check path.".format(self.input), display=True)
            sys.exit(1)
            
        if (self.snp_col != "" and self.snp_col not in df.columns) or (self.snp_col == ""):
            self.snp_col = df.columns[0]
            self.log("Warning: Defaulting to first column as SNP column.", display=True)
            
        if (self.allele_col != "" and self.allele_col not in df.columns) or (self.allele_col == ""):
            self.allele_col = df.columns[1]
            if df[self.allele_col].dtype != str:
                self.log("Error: Allele column is non-string. Please check the correct column is being read.", display=True)
                sys.exit(1)
            self.log("Warning: Defaulting to second column as allele column.", display=True)
            
        df['freq'] = -1
        
        loop = asyncio.get_event_loop()
        future = asyncio.ensure_future(self.get_data_asynchronous(df))
        loop.run_until_complete(future)
        
        df['freq'] = df['MarkerName'].map(self.d)
        df['freq'] = df['freq'].fillna(-1)
        
        try:
            df.to_csv(self.output, sep=self.osep, na_rep=None, index=False)
        except Exception as e:
            self.log("Error {}: Could not output file at \"{}\".".format(e, self.output), display=True)
            sys.exit(1)
            
    def fetch(self, session, snp, allele):
        """
        Processes the request and parses it.
        Results are added to the dictionary which is later merged into the
        input file.

        Parameters
        ----------
        session : Session object
            This is the session for the requests.
        snp : str
            Full SNP name. SNP must be of "rs" type and is clipped here.
        allele : str
            Allele whose frequency will be looked up.

        Returns
        -------
        None.

        """
        base_url = "https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/"
        for _ in range(self.attempts):
            try:
                with requests.get(base_url + snp[2:] + '/frequency', timeout=self.timeout) as response:
                    if response.status_code != 200:
                        self.d[snp] = -1
                        response.close()
                        return
                    
                    data = response.json()
                    response.close()
                    try:
                        res = data['results']
                    except KeyError:
                        self.d[snp] = -1
                        return
                    
                    try:
                        res = res[list(res.keys())[0]]['counts'][self.organism]['allele_counts'][self.population]
                    except KeyError:
                        self.log("Warning: SNP {} could not be found. Likely due to organism ({})/population ({}) not found.".format(snp, self.organism, self.population), display=False)
                        self.d[snp] = -1
                        return
        
                    alleles = list(res.keys())
                    counts = list(res.values())
                    sum_ = sum(counts)
                    if sum_ != 0:
                        for i in range(len(alleles)):
                            if alleles[i].upper() == allele.upper():
                                self.d[snp] = counts[i] / sum_
                                #print(snp, counts[i]/sum_)
                                break
                    return
            except requests.exceptions.Timeout:
                self.d[snp] = -1
                self.log("Timeout for {}, attempt: {}".format(snp, _), display=False)
                pass
            continue
            
    async def get_data_asynchronous(self, df):
        """
        Asynchronous function which handles the requests sessions.

        Parameters
        ----------
        df : DataFrame
            Dataframe of the input file which is cut to the SNPs and allele
            columns.

        Returns
        -------
        None.

        """
        snps = df[self.snp_col]
        alleles = df[self.allele_col]
        reference = dict(zip(snps, alleles))
        
        with ThreadPoolExecutor(max_workers = self.workers) as executor:
            with requests.Session() as session:
                loop = asyncio.get_event_loop()
                
                tasks = [loop.run_in_executor(executor, self.fetch, *(session, snp, allele)) for snp, allele in reference.items()]
                
                for response in await asyncio.gather(*tasks):
                    pass
            
    def log(self, msg, display=False):
        """
        Logs a message.

        Parameters
        ----------
        msg : str
            Message to be output to user.
        display : bool, optional
            Should this be printed as well as writen to the log file.
            The default is False.

        Returns
        -------
        None.

        """
        display = (self.verbose == 0 or display == True)
        if display:
            print(msg)
            
        try:
            f = open("eaf.log", 'a')
            f.write(msg+"\n")
            f.close()
        except:
            pass
"""        
def main(input_f, 
            isep="\t", 
            output_f="alff_output.txt", 
            osep="\t",
            snp_col="",
            allele_col="",
            organism="PRJNA507278",
            population="SAMN10492695",
            workers=5,
            timeout=5,
            attempts=5,
            verbose=False):
"""
"""
    This function will allow use of this module as an import.
    Wraps arguments into sys.argv and passes these to the ALFF class
    call and returns as normal.

    Parameters
    ----------
    input_f : str
        Path to input file.
    isep : str, optional
        Separator of input file. The default is "\t".
    output_f : str, optional
        Path to output file. The default is "alff_output.txt".
    osep : str, optional
        Separator of output file. The default is "\t".
    snp_col : str, optional
        Name of the SNP column in the input file. If no name is given, this
        will default to the first column.
    allele_col : str, optional
        Name of the allele column in the input file for which frequencies will
        be found. If no name is given, this will default to the second column.
    organism : str, optional
        Organism for whose allele frequencies to search. Please see 
        https://www.ncbi.nlm.nih.gov/bioproject/browse for more details. 
        The default is "PRJNA507278" (or human).
    population : str, optional
        Population subtype to use. Please see
        https://www.ncbi.nlm.nih.gov/snp/docs/gsr/alfa/ for more details. 
        The default is "SAMN10492695" which is European.
    workers : int, optional
        How many requests to send at once. Be mindful to not overload servers 
        by setting too high. Setting this lower will take longer but is less 
        likely to timeout. The default is 5.
    timeout : int, optional
        (Seconds) How long to wait before a request is considered as timed out.
        The default is 5.
    attempts : int, optional
        How many attempts should be tried if the request is met with a timeout 
        or other related error. The default is 5.
    verbose : bool, optional
        (Boolean) Will print verbose warnings and messages. This will be 
        automatically written to 'eaf.log' in the folder where this is run. 
        The default is False.

    Returns
    -------
    ALFF class.

"""
"""
    args = ['-i', input_f,
            '-is', isep,
            '-o', output_f,
            '-os', osep,
            '-snp', snp_col,
            '-allele', allele_col,
            '-org', organism,
            '-pop', population,
            '-w', workers,
            '-t', timeout,
            '-a', attempts,
            '-v', verbose]
"""

if __name__ == "__main__":
    # If this script is run in iPython there is an error that will cause the
    # asyncio loop to not run as it is run already. This fixes that.
    try:
        shell = get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell' or shell == 'TerminalInteractiveShell':
            nest_asyncio.apply()
    except NameError:
        pass
    
    _start_time = default_timer()
    
    app = ALFF()
    app.run()
    
    _time_taken = "{:5.2f}s".format(default_timer() - _start_time)
    app.log("Info: Time taken to run: {}".format(_time_taken), display=True)
    sys.exit(0)
    