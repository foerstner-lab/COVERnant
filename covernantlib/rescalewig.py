def rescale_wiggle(args):
    with open(args.output_prefix + ".wig", "w") as output_fh:
        with open(args.input_file) as input_fh:
            for line in input_fh:
                if not (line.startswith("variableStep")
                        or line.startswith("track")
                        or len(line.strip()) == 0):
                    pos, value = line.split()
                    line = "{} {}\n".format(
                        pos, float(value) * args.factor)
                output_fh.write(line)
