import sys
from datetime import datetime
import math
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def plot_matrix(args):
    pp = PdfPages(args.output_file)
    matrix = pd.read_table(args.matrix_file)
    if args.share_y_max:
        max_y = np.nanmax([
            np.nanmax(row[4:]) for row in map(list, matrix.values)])
    x_values = [int(x) for x in list(matrix.columns)[4:]]
    if args.share_x_range:
        min_x = np.nanmin(x_values)
        max_x = np.nanmax(x_values)
    for row in map(list, matrix.values):
        range_name = "{} {} {} {}".format(row[0], row[1], row[2], row[3])
        log("Plotting {}".format(range_name))
        y_values = row[4:]
        assert len(x_values) == len(y_values)
        filtered_x_values = []
        filtered_y_values = []
        for x, y in zip(x_values, y_values):
            if not math.isnan(y):
                filtered_x_values.append(x)
                filtered_y_values.append(y)
        if args.share_x_range:
            cur_min_x = min_x
            cur_max_x = max_x
        else:
            cur_min_x = np.nanmin(filtered_x_values)
            cur_max_x = np.nanmax(filtered_x_values)
        cur_max_y = max_y if args.share_y_max else np.nanmax(filtered_y_values)
        fig = plt.figure()
        plt.plot(x_values, y_values, color="black", )
        plt.title(range_name)
        plt.axis([cur_min_x, cur_max_x, 0, cur_max_y])
        plt.xlabel("Position")
        plt.ylabel("Coverage")
        pp.savefig()
        plt.close(fig)
    pp.close()

    
def log(msg):
    sys.stderr.write("{} - {}\n".format(
        datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S'), msg))
