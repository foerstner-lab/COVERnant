COVERnant
=========

COVERnant is a tool for the generation and manipulation of coverage
files (currently in wiggle format) of high-throughput sequencing data.

The tool is currently in an **early development stage**.

COVERnant has several subcommands as its command line help shows:

::

    $ covernant -h
    usage: covernant [-h] [--version]
                     {ratio,extract,plot_matrix,bed_to_wig,rescale_wig} ...

    positional arguments:
      {ratio,extract,plot_matrix,bed_to_wig,rescale_wig}
                            commands
        ratio               Generate ratio plots of two alignment files in Bam
                            formar.
        extract             Extract coverage values from a wiggle file based on
                            coordinates in a bed file and generate a matrix.
        plot_matrix         Plot the content of the extracted coverage matrix.
        bed_to_wig          Converts Bed files to coverage files in wiggle formats
        rescale_wig         Multiplies each value of a wiggle file with a given
                            factor.

    optional arguments:
      -h, --help            show this help message and exit
      --version, -v         show version

Subcommand ``ratio``
--------------------

::

    usage: covernant ratio [-h] [--denominator DENOMINATOR_BAM_FILE]
                           [--numerator NUMERATOR_BAM_FILE]
                           [--output_prefix OUTPUT_PREFIX] [--paired_end]
                           [--window_size WINDOW_SIZE] [--step_size STEP_SIZE]
                           [--factor_numerator FACTOR_NUMERATOR]
                           [--factor_denominator FACTOR_DENOMINATOR]
                           [--keep_zero_coverage]
                           [--denominator_name DENOMINATOR_NAME]
                           [--numerator_name NUMERATOR_NAME]
                           [--ratio_name RATIO_NAME]

    optional arguments:
      -h, --help            show this help message and exit
      --denominator DENOMINATOR_BAM_FILE
      --numerator NUMERATOR_BAM_FILE
      --output_prefix OUTPUT_PREFIX, -o OUTPUT_PREFIX
      --paired_end          Paired reads are treated as one fragment an the start
                            and end positions are used accordingly
      --window_size WINDOW_SIZE
                            Window size for sliding window average calculation.
                            Must be an odd number. (Default is 1).
      --step_size STEP_SIZE
                            Step size for sliding window average calculation.
                            Default is 1.
      --factor_numerator FACTOR_NUMERATOR
                            A factor the numerator values are are multiplied with.
      --factor_denominator FACTOR_DENOMINATOR
                            A factor the denominator values are are multiplied
                            with.
      --keep_zero_coverage  Also write coordinates that have a coverage of 0.
                            Default is to discard those.
      --denominator_name DENOMINATOR_NAME
      --numerator_name NUMERATOR_NAME
      --ratio_name RATIO_NAME

Subcommand ``extract``
----------------------

::

    $ covernant extract -h
    usage: covernant extract [-h] [--output_prefix OUTPUT_PREFIX]
                             [--flip_reverse_strand]
                             [--matrix_alignment {left,center,right}]
                             [--window_size WINDOW_SIZE] [--step_size STEP_SIZE]
                             coverage_file coordinate_file

    positional arguments:
      coverage_file
      coordinate_file

    optional arguments:
      -h, --help            show this help message and exit
      --output_prefix OUTPUT_PREFIX
      --flip_reverse_strand
                            Flip the coverage value list of entries located at the
                            minus strand
      --matrix_alignment {left,center,right}
                            default is 'left'.
      --window_size WINDOW_SIZE
                            Window size for sliding window average calculation.
                            Must be an odd number.
      --step_size STEP_SIZE
                            Step size for sliding window average calculation.
                            Default is 1.

Subcommand ``plot_matrix``
--------------------------

::

    $ covernant plot_matrix -h
    usage: covernant plot_matrix [-h] [--output_file OUTPUT_FILE]
                                 [--share_x_range] [--share_y_max]
                                 matrix_file

    positional arguments:
      matrix_file

    optional arguments:
      -h, --help            show this help message and exit
      --output_file OUTPUT_FILE
      --share_x_range       Use the same x range in all plots.
      --share_y_max         Use the same maximum y value in all plots.

Subcommand ``bed_to_wig``
-------------------------

::

    $ covernant bed_to_wig -h
    usage: covernant bed_to_wig [-h] [--output_prefix OUTPUT_PREFIX]
                                [--window_size WINDOW_SIZE]
