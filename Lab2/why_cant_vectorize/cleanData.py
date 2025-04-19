def filter_lines_with_spread_functions(input_file, output_file=None):
    if output_file is None:
        output_file = input_file  # Overwrite the original file

    with open(input_file, 'r') as f:
        lines = f.readlines()

    filtered_lines = [line for line in lines if 'spread_functions' in line]

    with open(output_file, 'w') as f:
        f.writelines(filtered_lines)


filter_lines_with_spread_functions('not_vectorized_srcVect_-O1.log')
filter_lines_with_spread_functions('not_vectorized_srcVect_-O2.log')
filter_lines_with_spread_functions('not_vectorized_srcVect_-O3.log')
filter_lines_with_spread_functions('not_vectorized_srcVect_-Ofast.log')
