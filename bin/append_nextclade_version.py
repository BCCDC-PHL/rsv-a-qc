import argparse
import csv

def add_nextclade_version(input_file, output_file, nextclade_version):

    with open(input_file, 'r', newline='') as f_in:

        reader = csv.reader(f_in, delimiter='\t')
    
        header = next(reader)
    
        header.append('nextclade_version')
        
        with open(output_file, 'w', newline='') as f_out:

            writer = csv.writer(f_out, delimiter='\t')

            writer.writerow(header)
            
            for row in reader:
 
                row.append(nextclade_version)
   
                writer.writerow(row)

def main(args):
    add_nextclade_version(args.nextclade_tsv, args.output_name, args.nextclade_version)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--nextclade_tsv")
    parser.add_argument("--output_name")
    parser.add_argument("--nextclade_version")
    args = parser.parse_args()
    main(args)
