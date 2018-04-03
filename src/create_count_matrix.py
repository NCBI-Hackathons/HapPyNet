import argparse
import glob
import logging
import pandas as pd
import re
import sys

COUNT_FILE_EXTENSION = '.count'
COUNT_FILE_COLUMNS = ['chr', 'start', 'stop', 'count']
# Used to recover filename from input path.
# Must extend if you use input other than SRR.
EXPECTED_FILE_PREFIXES = ['SRR']


def get_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        '-d', '--debug',
        help='Print debug messages.',
        action='store_const',
        dest='log_level',
        const=logging.DEBUG,
        default=logging.WARNING
    )
    parser.add_argument(
        '-i', '--in_dir',
        help='The folder containing the .count files.',
        default='../results/counts/',
    )
    parser.add_argument(
        '-t', '--output_tag',
        help='Tag appended to output count matrix file.',
        default=''
    )
    return parser.parse_args()


def get_files_from_folder(folder):
    '''
    Searches folder for .count files.
    '''
    logging.debug('Checking {} for {} files'.format(folder, COUNT_FILE_EXTENSION))
    return glob.glob('{}*{}'.format(folder, COUNT_FILE_EXTENSION))


def load_df(f):
    '''
    Expects .count TSV sorted by chr column.
    '''
    return pd.read_table(f, header=None, names=COUNT_FILE_COLUMNS)


def parse_filename(filename):
    p = re.compile(
        '([{}]\w+)\{}'.format(
            '|'.join(EXPECTED_FILE_PREFIXES),
            COUNT_FILE_EXTENSION
        )
    )
    try:
        # Returns just the filename starting with e.g. SRR.
        return p.search(filename).group(1)
    except AttributeError:
        logging.warn(
            'Unable to recover any {} filenames from {}. Exiting'.format(
                '/'.join(EXPECTED_FILE_PREFIXES),
                filename
            )
        )
        sys.exit(1)


def load_dfs(files):
    '''
    Returns list of filename, dataframe tuples.
    '''
    return [(parse_filename(f), load_df(f)) for f in files]


def rename_column(df, old, new, append='_count'):
    '''
    Associate count column with sample name.
    '''
    return df.rename(columns={old: new + append})


def combine_dfs(dfs):
    '''
    Renames column by filename and concatenates count columns from dfs.
    '''
    dfs = [rename_column(d[1], 'count', d[0]) for d in dfs]
    dfs = pd.concat(dfs, axis=1)
    return rename_column(
        dfs[[c for c in dfs.columns if 'count' in c]].reset_index(),
        'index',
        'ld_block',
        append=''
    )


def output_matrix(matrix, tag):
    matrix.to_csv('haplo_phen_count_matrix{}.tsv'.format(tag), sep='\t', index=False)


def main():
    '''
    Returns count matrix of LD block by sample given folder of .count files.'
    '''
    args = get_args()
    logging.basicConfig(level=args.log_level)
    count_files = get_files_from_folder(args.in_dir)
    if count_files:
        logging.debug('Found {} files'.format(count_files))
        dfs = load_dfs(count_files)
    else:
        logging.warn('No count files detected. Exiting.')
        sys.exit()
    matrix = combine_dfs(dfs)
    output_matrix(matrix, args.output_tag)


if __name__ == '__main__':
    main()
