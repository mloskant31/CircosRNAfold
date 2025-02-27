from argparse import ArgumentParser
import pyBigWig
import pybedtools
import pycircos
import RNA
import collections
import cv2
import os
import math


# Command line flags
parser = ArgumentParser()
parser.add_argument("-g", "--genes", dest="genes", help="The genomic coordinates of genes in BED format", required=True)
parser.add_argument("-b", "--bindings", dest="bindings", nargs='+',
                    help="The binding sites of proteins in BED format (max. 4 files)",
                    required=False)
parser.add_argument("-fi", "--fasta", dest="fasta", help="The chromosomal sequences as a FASTA file", required=True)
parser.add_argument("-o", "--output", dest="output", help="Name of the folder in which all PNGs are saved", required=False,
                    default="myCircos")
parser.add_argument("-bw", dest="bw", nargs='+', help="BigWig file(s) with the raw data (max. 4 files)", required=False)
parser.add_argument("-color", dest="color", nargs='+',
                    help="The colors for the iCLIP plots. If none or too few/too many colors are specified " +
                         "than there are specified bw files, default colors are used (max. 4 colors)", required=False)
parser.add_argument("-m", dest="mature", help="The mature miRNA locations in BED format", required=False)


args = parser.parse_args()
genes = str(args.genes)
fasta_file = args.fasta
folder = str(args.output)
mirna_mature = str(args.mature)

# create directory for circos plots
os.mkdir(folder)

if args.bindings and len(args.bindings) > 4:
    parser.error("A maximum of 4 arguments may be specified for --bindings.")
if len(args.bw) > 4:
    parser.error("A maximum of 4 arguments may be specified for -bw.")

if args.color is None or len(args.color) != len(args.bw):
    print("Not as many colors as bw files were specified or no colors, therefore default colors are used")
    color_bws = ["#800080", "#006400", "#FF6347", "#4682B4"]  # purple, green, red, blue
    color_bws = color_bws[:len(args.bw)]
else:
    color_bws = args.color
    color_bws = ["#" + color for color in color_bws]


def convert_dot_bracket(structure, gene_name):
    """
    Converts dot-bracket structure for circos

    Parameters
    ----------
    structure : str
        Dot-bracket structure of a sequence
    gene_name: str
        Name of the gene

    Returns
    -------
    start_list: list
        Start data location of linked data
    end_list: list
        End data location of linked data
    """
    start_stack = []  # Stack to track open brackets
    start_list = []
    end_list = []

    for position, char in enumerate(structure):
        if char == "(":
            start_stack.append(position)
        elif char == ")":
            if len(start_stack) > 0:
                start = start_stack.pop()  # Pop the last open bracket position
                start_list.append((gene_name, start, start, 499))
                end_list.append((gene_name, position, position, 499))

    return start_list, end_list


def find_features_on_gene(bed_file, gene_name, gene_chrom, gene_chromstart, gene_chromend):
    """
    Finds and saves the start positions and widths of the features which are on the current gene

    Parameters
    ----------
    bed_file: str
        BED file
    gene_name: str
        Name of the gene
    gene_chrom: str
        Chromosome on which the gene is located
    gene_chromstart: int
        Start of the sequence of the gene
    gene_chromend: int
        End of the sequence of the gene

    Returns
    -------
    arcdata_dict: dict
        A dictionary of the start positions and widths of the sequence on the circo plot
    """
    arcdata_dict = collections.defaultdict(dict)
    bed_file = pybedtools.BedTool(bed_file)
    # Filter out all irrelevant binding sites
    for feature in bed_file:
        if feature.chrom == gene_chrom:
            if feature.start <= gene_chromend and feature.end >= gene_chromstart:
                if feature.strand == gene_strand:
                    if feature.strand == "-":
                        start = int(feature.end) - gene_chromend
                    else:
                        start = int(feature.start) - gene_chromstart
                    width = int(feature.end) - int(feature.start)
                    if gene_name not in arcdata_dict:
                        arcdata_dict[gene_name]["positions"] = []
                        arcdata_dict[gene_name]["widths"] = []
                    # Save position on the gene
                    arcdata_dict[gene_name]["positions"].append(abs(start))
                    # Save length of the binding site
                    arcdata_dict[gene_name]["widths"].append(width)
    return arcdata_dict


def hex_to_rgb(hex_color):
    """
    Converts a hexadecimal color string to an RGB tuple.

    Parameters
    ----------
    hex_color: str
        A string representing a color in hexadecimal format.

    Returns
    -------
    rgb: tuple
        A tuple containing three integers that represent the red, green, and blue components of the color.
    """
    hex_color = hex_color.lstrip('#')
    rgb = tuple(int(hex_color[i:i + 2], 16) for i in (0, 2, 4))
    return rgb


# Load the BED file
bed_file = pybedtools.BedTool(genes)

# Run through the entries in the BED file
for feature in bed_file:
    gene_sequence = ""
    gene_chrom = feature.chrom  # chromosome on which the gene is located
    gene_chromstart = feature.start  # start coordinate on the chromosome
    gene_chromend = feature.end  # end coordinate on the chromosome
    gene_strand = feature.strand  # strand orientation
    gene_id = feature.name  # id
    # Extract the sequence from the FASTA file
    sequence = bed_file.sequence(fi=fasta_file, name=True, s=True)
    with open(sequence.seqfn, 'r') as file:
        for line in file:
            # save the sequence at the appropriate line until the next gene begins
            if line.startswith(">" + gene_id):
                try:
                    line = next(file)
                    gene_sequence = gene_sequence + line
                    line = next(file)
                    while not (line.startswith(">")):
                        gene_sequence = gene_sequence + line
                        line = next(file)
                except StopIteration:
                    break

    # Folding of the RNA sequence
    # Package viennarna version: 2.6.4, RNAfold
    sequence = gene_sequence.strip()
    structure, _ = RNA.fold(sequence)

    # Circos plot
    Garc = pycircos.Garc
    Gcircle = pycircos.Gcircle
    circle = Gcircle(figsize=(25, 25))
    # The circle itself
    label = gene_id
    arc = Garc(arc_id=gene_id, size=len(sequence), raxis_range=(500, 500), labelposition=460, labelsize=30,
               label_visible=True, edgecolor="black", facecolor="black", linewidth=1.5, label=label)
    circle.add_garc(arc)
    # The arc rectangle with start and end angle of the circos plot
    circle.set_garcs(0, 360)

    # Markings on the circle
    # BED file starts with 0, so does the circos plot
    if len(sequence) > 2500:
        tick_labels = list(range(0, len(sequence), 1000))
    else:
        tick_labels = list(range(0, len(sequence), 100))
    circle.tickplot(gene_id, raxis_range=(500, 510), tickwidth=2.0,
                    tickpositions=tick_labels, ticklabels=tick_labels,
                    ticklabeldirection="outer")

    # Tick plot for sequence
    # Only printed if the sequence is not too long for the bases to be readable on the plot
    # No longer visible if there are more than two bw files
    if len(sequence) <= 240 and len(args.bw) < 3:
        rna_sequence = sequence.replace("T", "U")
        tick_positions = list(range(0, len(rna_sequence), 1))
        circle.tickplot(gene_id, raxis_range=(790, 80), tickwidth=0.0,
                        tickpositions=tick_positions,
                        ticklabels=rna_sequence,
                        ticklabeldirection="inner")

    # Chord plot for each base pairing
    # Converting dot-bracket structure for circos
    start_list, end_list = convert_dot_bracket(structure, gene_id)
    for i in range(0, len(start_list)):
        circle.chord_plot(start_list[i], end_list[i], linewidth=1, edgecolor="grey")

    # Barplot for raw data
    raxis_range_start = 603
    raxis_range_end = raxis_range_start + 60
    if args.bw is not None:
        bw = args.bw
        for bw_index, color_bws_index in zip(bw, color_bws):
            current_bw = pyBigWig.open(bw_index)
            # Converting every nan data to zero data
            values_without_nan = [0 if math.isnan(x) else x for x in
                                  current_bw.values(gene_chrom, gene_chromstart, gene_chromend)]
            if gene_strand == "-":
                # mirrored value for minus strand gene
                values_without_nan = values_without_nan[::-1]
            # If the data set consists only of zeros, no plot is created
            if not all(value == 0 for value in values_without_nan):
                circle.barplot(gene_id, data=values_without_nan,
                               positions=list(range(0, len(sequence), 1)),
                               width=1, raxis_range=[raxis_range_start, raxis_range_end],
                               edgecolor=color_bws_index, facecolor=color_bws_index, linewidth=0.5, spine=False)

            raxis_range_start = raxis_range_end + 22
            raxis_range_end = raxis_range_start + 60

    # Barplot for protein bindings
    if args.bindings is not None:
        # Position of the bindings
        raxis_range_start = 580
        raxis_range_end = raxis_range_start + 20

        bindings = args.bindings
        color = ["#FF69B4", "#ADFF2F", '#FFB3BA', '#B3D9FF']
        color_index = 0

        for binding in bindings:

            arcdata_dict = find_features_on_gene(binding, gene_id, gene_chrom, gene_chromstart, gene_chromend)
            for key in arcdata_dict:
                circle.barplot(key, data=[1] * len(arcdata_dict[key]["positions"]),
                               positions=arcdata_dict[key]["positions"],
                               width=arcdata_dict[key]["widths"], raxis_range=[raxis_range_start, raxis_range_end],
                               edgecolor="black", facecolor=color[color_index], linewidth=0.5)

            color_index += 1  # next color
            raxis_range_start = raxis_range_end + 62
            raxis_range_end = raxis_range_start + 20

            # Barplot for mature miRNA
        if args.mature is not None:
            mirna_color = "#FF8800"
            mirna_mature = args.mature
            mature_dict = find_features_on_gene(mirna_mature, gene_id, gene_chrom, gene_chromstart, gene_chromend)
            for key in mature_dict:
                circle.barplot(key, data=[1] * len(mature_dict[key]["positions"]),
                               positions=mature_dict[key]["positions"],
                               width=mature_dict[key]["widths"],
                               raxis_range=[raxis_range_start, raxis_range_end],
                               edgecolor="black",
                               facecolor=mirna_color)

        # Save the plot as specified in the command line
        circle.figure.savefig(folder + "/" + gene_id + ".png")

        # Legende
        # Position of the legend/rectangle
        position_y = 10
        position_y2 = 40
        position_x = 10

        # Load the image
        image = cv2.imread(folder + "/" + gene_id + ".png")

        # Font, size, color and thickness for the text
        font = cv2.FONT_HERSHEY_SIMPLEX
        font_scale = 1.1
        font_color = (0, 0, 0)
        font_thickness = 2
        legend_position = (position_x, position_y2)
        if args.mature is not None:
            # rectangle (BGR color)
            rgb = hex_to_rgb(mirna_color)
            rectangle_color = (rgb[2], rgb[1], rgb[0])
            rectangle_top_left = (position_x, position_y)
            rectangle_bottom_right = (position_x + 40, position_y2)
            cv2.rectangle(image, rectangle_top_left, rectangle_bottom_right, rectangle_color, thickness=cv2.FILLED)

            legend_position = (position_x + 45, position_y2)
            cv2.putText(image, "miRNA(*)", legend_position, font, font_scale, font_color, font_thickness)

        position_y += 50
        position_y2 += 50

        legend_position = (position_x, position_y2)
        cv2.putText(image, "Binding sites from", legend_position, font, font_scale, font_color, font_thickness)

        position_y += 50
        position_y2 += 50
        color_index = 0

        for binding in bindings:
            legend_text = binding[:-4]

            # rectangle (BGR color)
            rgb = hex_to_rgb(color[color_index])
            rectangle_color = (rgb[2], rgb[1], rgb[0])
            rectangle_top_left = (position_x, position_y)
            rectangle_bottom_right = (position_x + 40, position_y2)
            cv2.rectangle(image, rectangle_top_left, rectangle_bottom_right, rectangle_color, thickness=cv2.FILLED)

            legend_position = (position_x + 45, position_y2)
            cv2.putText(image, legend_text, legend_position, font, font_scale, font_color, font_thickness)

            # position for next text and rectangle
            position_y += 40
            position_y2 += 40
            # next color
            color_index += 1

        position_y += 5
        position_y2 += 5

        legend_position = (position_x, position_y2)
        cv2.putText(image, "Crosslink data for", legend_position, font, font_scale, font_color, font_thickness)

        position_y += 50
        position_y2 += 50

        for i, bw_index in zip(color_bws, bw):
            rgb = hex_to_rgb(i)
            rectangle_color = (rgb[2], rgb[1], rgb[0])
            rectangle_top_left = (position_x, position_y)
            rectangle_bottom_right = (position_x + 40, position_y2)
            cv2.rectangle(image, rectangle_top_left, rectangle_bottom_right, rectangle_color, thickness=cv2.FILLED)

            # Font, size, colour and thickness for the text
            legend_position = (position_x + 45, position_y2)
            cv2.putText(image, bw_index[:-3], legend_position, font, font_scale, font_color, font_thickness)
            position_y += 40
            position_y2 += 40

        # Save the image with the added legend
        cv2.imwrite(folder + "/" + gene_id + ".png", image)

