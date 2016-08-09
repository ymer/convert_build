import pandas as pd
import argparse
import subprocess
import os


def parse_input():
    args = {'data': 'plink files (filestub)',
            'posonref': 'NCBI SNPChrPosOnRef file',
            'out': 'output filenames (filestub)'}

    parser = argparse.ArgumentParser()
    for arg, hlp in args.items():
        parser.add_argument('--' + arg, action='store', dest=arg, help=hlp)

    return parser.parse_args()


s = parse_input()

data = pd.read_csv(s.data + '.bim', header=None, sep='\s+', names=['chrom_orig', 'snp', 'gpos', 'pos_orig', 'A1', 'A2'])
print('Data read. {} snps read.'.format(len(data)))

posonref = pd.read_csv(s.posonref, header=None, sep='\s+', names=['snp', 'chrom', 'pos', 'strand'])
posonref.snp = posonref.snp.apply(lambda x: 'rs' + str(int(x)))
posonref.chrom = posonref.chrom.astype(str)
posonref.chrom.replace({'AltOnly': 25, 'PAR': 25, 'NotOn': 25, 'Un': 25, 'Y': 24, 'X': 23, 'MT': 26}, inplace=True)
print('Posonref read. {} snp positions read.'.format(len(posonref)))

data = pd.merge(data, posonref, on='snp', how='left')
data = data[['chrom', 'snp', 'gpos', 'pos', 'A1', 'A2']]
print('Positions replaced with values from posonref.')

keepsnps = data.dropna()[['snp']]
keepsnps.to_csv('keepsnps.txt', header=False, index=False)
print('{} snps not found inf posonref. These will be removed removed'.format(len(data) - len(keepsnps)))

data.fillna(0, inplace=True)
data[['chrom', 'pos']] = data[['chrom', 'pos']].astype(int, inplace=True)
data.to_csv('temp.bim', header=False, index=False, sep=' ')

subprocess.call(['plink', '--bim', 'temp.bim', '--fam', s.data + '.fam', '--bed', s.data + '.bed', '--extract', 'keepsnps.txt', '--make-bed', '--out', s.out])

os.remove('temp.bim')
os.remove('keepsnps.txt')

