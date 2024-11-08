import csv
import argparse


def convert_ucsc_genome_bed_to_refseq(ucsc_genome_bed_file, sequence_report_file, out_file):

    conversion_dict = {}

    with open(sequence_report_file, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            conversion_dict[row[11]] = row[8]

    with open(ucsc_genome_bed_file, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        with open(out_file, "w") as out:
            writer = csv.writer(out, delimiter="\t")
            for row in reader:
                if row[0] in conversion_dict:
                    row[0] = conversion_dict[row[0]]
                    writer.writerow(row)
                else:
                    print(f"Error: {row[0]} not found in conversion dictionary")
                    continue

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--ucsc_genome_bed_file", help="Path to UCSC genome bed file")
    parser.add_argument("-r", "--sequence_report_file", help="Path to sequence report file")
    parser.add_argument("-o", "--out_file", help="Path to output file")
    args = parser.parse_args()

    convert_ucsc_genome_bed_to_refseq(args.ucsc_genome_bed_file, args.sequence_report_file, args.out_file)
