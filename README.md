# COVERnant

A tool to generate and manipulate coverage plots of high-throughput sequencing data.


```
usage: covernant [-h] [--version] {ratio,extract,plot_matrix,bed_to_wig} ...

positional arguments:
  {ratio,extract,plot_matrix,bed_to_wig}
                        commands
    ratio               Generate ratio plots of two alignment files.
    extract             Extract coverage values from wiggle file.
    plot_matrix         Plot the content of the extracted coverage matrix.
    bed_to_wig          Converts Bed file to coverage in wiggle formats

optional arguments:
  -h, --help            show this help message and exit
  --version, -v         show version
```


## Subcommand `ratio`

```
$ covernant ratio -h
usage: covernant ratio [-h] [--output OUTPUT_PREFIX] [--paired_end]
                       [--window_size WINDOW_SIZE] [--step_size STEP_SIZE]
                       [--factor FACTOR] [--keep_zero_coverage]
                       [--denominator_name DENOMINATOR_NAME]
                       [--numerator_name NUMERATOR_NAME]
                       [--ratio_name RATIO_NAME]
                       denominator_bam_file numerator_bam_file
covernant ratio: error: the following arguments are required: denominator_bam_file, numerator_bam_file


covernant ratio -h
usage: covernant ratio [-h] [--output OUTPUT_PREFIX] [--paired_end]
                       [--window_size WINDOW_SIZE] [--step_size STEP_SIZE]
                       [--factor FACTOR] [--keep_zero_coverage]
                       [--denominator_name DENOMINATOR_NAME]
                       [--numerator_name NUMERATOR_NAME]
                       [--ratio_name RATIO_NAME]
                       denominator_bam_file numerator_bam_file

positional arguments:
  denominator_bam_file
  numerator_bam_file

optional arguments:
  -h, --help            show this help message and exit
  --output OUTPUT_PREFIX, -o OUTPUT_PREFIX
  --paired_end          Paired reads are treated as one fragment an the start
                        and end positions are used accordingly
  --window_size WINDOW_SIZE
                        Window size for sliding window average calculation.
                        Must be an odd number. (Default is 1).
  --step_size STEP_SIZE
                        Step size for sliding window average calculation.
                        Default is 1.
  --factor FACTOR       A factor the final ratio is multiplied with.
  --keep_zero_coverage  Also write coordinates that have a coverage of 0.
                        Default is to discard those.
  --denominator_name DENOMINATOR_NAME
  --numerator_name NUMERATOR_NAME
  --ratio_name RATIO_NAME
```

## Subcommand `extract`

```
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
```

## Subcommand `plot_matrix`

```
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
```

## Subcommand `bed_to_wig`

```
$ covernant bed_to_wig -h
usage: covernant bed_to_wig [-h] [--output_prefix OUTPUT_PREFIX]
                            [--window_size WINDOW_SIZE]

```


