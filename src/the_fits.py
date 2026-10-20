

import argparse

import pandas as pd


REQUIRED_COLUMNS = {
    'c_lowest_prob',
    'in_c_lowest_cluster',
    'c1_max',
    'tumor_VAF',
    # lowest_correction groups rows by sample when updating c1_max.
    'Sample',
}


def lowest_correction(df):
    '''
    Corrects the lowest cluster probability and related fields for rows where c_lowest_prob is NaN and tumor_VAF is >= 1.0.
    
    Parameters:
    df (pd.DataFrame): The input DataFrame containing the columns 'c_lowest_prob', 'in_c_lowest_cluster', 'c1_max', and 'tumor_VAF'.

    Returns:
    pd.DataFrame: The corrected DataFrame.
    '''
    df.loc[(df.c_lowest_prob.isna()) & (df.tumor_VAF >= 1.0),'c_lowest_prob'] = 0.0
    df.loc[(df.c_lowest_prob.isna()) & (df.tumor_VAF >= 1.0),'in_c_lowest_cluster'] = False
    for sampleId in df.Sample.unique():
        c1_max = df[df.Sample == sampleId].c1_max.max()
        df.loc[(df.c_lowest_prob.isna()) 
                & (df.tumor_VAF >= 1.0) 
                & (df.Sample == sampleId),'c1_max'] = c1_max 
    return df


def correct_csv(input_csv, output_csv=None):
    """Correct a CSV and write it to the input path or an optional output path."""
    df = pd.read_csv(input_csv)
    missing_columns = REQUIRED_COLUMNS.difference(df.columns)
    if missing_columns:
        missing = ', '.join(sorted(missing_columns))
        raise ValueError(f"Input CSV is missing required columns: {missing}")

    df = lowest_correction(df)
    df.to_csv(output_csv or input_csv, index=False)
    return df


def main():
    parser = argparse.ArgumentParser(
        description='Correct lowest-cluster values in a CSV file.'
    )
    parser.add_argument('input_csv', help='CSV file to correct')
    parser.add_argument(
        '-o', '--output-csv',
        help='Path for the corrected CSV (defaults to overwriting the input file)',
    )
    args = parser.parse_args()
    correct_csv(args.input_csv, args.output_csv)


if __name__ == '__main__':
    main()
