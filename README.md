# ALFF (Allele Frequency Finder)
This simple Python script uses NCBI's [variation service](https://api.ncbi.nlm.nih.gov/variation/v0/) (in particular the ALFA project) to find allele frequencies for a given list of SNPs and alleles. The tool is designed for use from the command line.

## Getting Started

### Prerequisites

This tool requires Python 3 to be installed. A list of dependencies is given in `requirements.txt`.

### Installation
Simple download or clone this repo like you would normally. No further installation required. 

### Usage
```
ALFF.py -h
```
lists a detailed breakdown of arguments which may be used. This tool requires only a single argument for the input file, though it is recommended you observe the rest of the available options to tailor the tool to your data.

## License

This project is licensed under the GNU v3.0 License - see the [LICENSE.md](LICENSE.md) file for details.
