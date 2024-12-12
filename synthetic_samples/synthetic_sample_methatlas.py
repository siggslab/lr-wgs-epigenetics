import os.path
import pysam
import yaml
import sys

## Remove any samples with weighting 0, so its easier
def filter_samples(config):
    samples = config['subsampling']['samples']
    new_samples = {}
    for s in samples:
        sample = samples[s]
        if float(sample['weight']) > 0:
            new_samples[s] = samples[s]
    config['subsampling']['samples'] = new_samples
    return config

########################
## Subsample bamfiles ##
########################
def subsample_bamfiles(config):
    # Get read count per bamfile
    subsampling_per_bamfile_n_reads_in_each_bamfile = {}
    samples = config['subsampling']['samples']
    for s in config['subsampling']['samples']:
        sample = samples[s]
        n_reads = int(sample['n_reads'])
        if n_reads < 0:
            print(f"finding n_reads in {sample['name']}")
            # Then we need to actually check the read count. -1 is a placeholder for idk
            n_reads = 0
            bamfile = pysam.AlignmentFile(sample['path'], 'rb')
            for _ in bamfile.fetch():
                n_reads += 1
            bamfile.close()
            print(f"found {n_reads} reads")

        subsampling_per_bamfile_n_reads_in_each_bamfile[s] = n_reads

    # Get subsampling fraction for each bamfile
    print("calculating subsample fractions")
    bamfile_subsample_fractions = {}
    total_weight = sum([float(samples[s]['weight']) for s in samples])
    for s in samples:
        sample = samples[s]
        weight = float(sample['weight'])
        fraction_of_total_reads = weight / total_weight

        n_reads = int(fraction_of_total_reads * int(config['subsampling']['total_reads']))
        print(f"{s}: {weight}/{total_weight} -> {fraction_of_total_reads * 100:.5f}% {n_reads}/{config['subsampling']['total_reads']} reads")


        subsample_fraction = n_reads / subsampling_per_bamfile_n_reads_in_each_bamfile[s]
        bamfile_subsample_fractions[s] = subsample_fraction

    # Actually subsample each bamfile
    subsampling_directory = os.path.join(config['working_directory'], 'subsampling')
    os.makedirs(subsampling_directory, exist_ok=True)
    for s in samples:
        sample = samples[s]
        if bamfile_subsample_fractions[s] == 0:
            continue
        print(f"subsampling {sample['name']} for {int(bamfile_subsample_fractions[s] * float(sample['n_reads']))}/{sample['n_reads']} reads -> {bamfile_subsample_fractions[s]}")
        sample_bamname = sample['name'].rstrip('.bam') + '.bam'  # Name of the output file
        out_path = os.path.join(subsampling_directory, sample_bamname)  # Output file path

        # subsample
        # pysam.view('-b',
        #            '-o', out_path,
        #            '-s', f"{bamfile_subsample_fractions[s]}",
        #            '-@', f"{config['threads']}",
        #            sample['path'])
        ## ^ apparently the pysam.view function is broken and cant write a new file
        os.system(f"samtools view -b -o {out_path} -s {bamfile_subsample_fractions[s]} -@ {config['threads']} {sample['path']}")

    # Concatenate into one bamfile
    print("concatenating")
    # pysam.cat(os.path.join(subsampling_directory, '*'), '-o', os.path.join(subsampling_directory, 'subsampled.bam'))
    os.system(f"samtools cat {os.path.join(subsampling_directory, '*')} -o {os.path.join(subsampling_directory, 'subsampled.bam')}")

    # Sort
    print("sorting")
    pysam.sort(
        '-@', f"{config['threads']}",
        '-o', os.path.join(subsampling_directory, 'subsampled.sorted.bam'),
        os.path.join(subsampling_directory, 'subsampled.bam'))
    # Index
    print("indexing")
    pysam.index(os.path.join(subsampling_directory, 'subsampled.sorted.bam'))

#####################
## Bam to Pat/Beta ##
#####################
def bam_to_pat(config):
    print("Converting to pat/beta")
    subsampling_directory = os.path.join(config['working_directory'], 'subsampling')
    wgbstools_directory = os.path.join(config['working_directory'], 'wgbstools')
    os.makedirs(wgbstools_directory, exist_ok=True)
    os.system(f"wgbstools bam2pat -@ {config['threads']} -np -o {wgbstools_directory} {os.path.join(subsampling_directory, 'subsampled.sorted.bam')}")
    # os.system(f"mv {os.path.join(subsampling_directory, '*.pat.gz')} {wgbstools_directory}")
    # os.system(f"mv {os.path.join(subsampling_directory, '*.beta')} {wgbstools_directory}")

####################
## Beta to Blocks ##
####################
def beta_to_blocks(config):
    print("getting blocks")
    wgbstools_directory = os.path.join(config['working_directory'], 'wgbstools')
    os.system(f"wgbstools beta_to_blocks --blocks_file {os.path.join(wgbstools_directory, 'subsample.blocks')} -@ {config['threads']} {os.path.join(wgbstools_directory, 'subsampled.sorted.beta')}")

#############################
## Beta/Blocks to Coverage ##
#############################
def beta_to_table(config):
    print("Getting methylation fractions")
    wgbstools_directory = os.path.join(config['working_directory'], 'wgbstools')
    os.system(f"wgbstools beta_to_table {os.path.join(wgbstools_directory, 'subsample.blocks')} -@ {config['threads']} --digits 5 --betas {os.path.join(wgbstools_directory, 'subsampled.sorted.beta')} -o {os.path.join(wgbstools_directory, 'subsampled_coverage.tsv')}")

#######################################
## CpG Coverage to Illumina Coverage ##
#######################################
def coverage_to_illumina(config):
    print("converting to Illumina sites")
    wgbstools_directory = os.path.join(config['working_directory'], 'wgbstools')
    methatlas_directory = os.path.join(config['working_directory'], 'methatlas')
    os.makedirs(methatlas_directory, exist_ok=True)

    mapping_file = config['illumina_mapping_file']
    """
    cg08169020,chr14:68790171-68790172 - 2bp, 1CpGs: 20366646-20366647
    """
    coverage_file = os.path.join(wgbstools_directory, 'subsampled_coverage.tsv')
    """
    chr	start	end	startCpG	endCpG	SAMPLE1	SAMPLE2	
    chr1	10469	10472	1	3	0.77064	0.90741
    """

    atlas_file = os.path.join(methatlas_directory, "subsampled_illumina_fractions.csv")

    cpg_map = open(mapping_file, 'r')
    coverage = open(coverage_file, 'r')
    atlas = open(atlas_file, 'w')

    coverage_dict = {}
    header = coverage.readline()
    samples = header.split("\t")[5:]

    line = coverage.readline()
    while line != "":
        split = line.split("\t")
        cpg_start = int(split[3])
        cpg_end = int(split[4])
        for i in range(cpg_start, cpg_end):
            c = {}
            for j, sample in enumerate(samples):
                c[sample] = split[5 + j]
            coverage_dict[i] = c

        line = coverage.readline()

    line = cpg_map.readline()
    illumina_converage = {}
    while line != "":
        split = line.split(",")
        illumina = split[0]
        cpg = int(split[2].split(" ")[-1].split("-")[0])

        try:
            illumina_converage[illumina] = coverage_dict[cpg]
        except:
            print(f"no coverage for {illumina}. Continuing")
        line = cpg_map.readline()

    header = "CpGs," + ",".join(samples) + "\n"
    atlas.write(header)
    for illumina in illumina_converage:
        outline = illumina
        for sample in samples:
            outline += "," + illumina_converage[illumina][sample]
        # outline += "\n"
        atlas.write(outline)

####################################
## Illumina Coverage to methatlas ##
####################################
def illumina_to_methatlas(config):
    methatlas_directory = os.path.join(config['working_directory'], 'methatlas')
    os.system(f"python {os.path.join(config['methatlas_path'], 'deconvolve.py')} -a {os.path.join(config['methatlas_path'], 'reference_atlas.csv')} -r -o {methatlas_directory} {os.path.join(methatlas_directory, 'subsampled_illumina_fractions.csv')}")


if __name__ == '__main__':
    with open(sys.argv[1], 'r') as f:
        config = yaml.safe_load(f)
    config = filter_samples(config)
    subsample_bamfiles(config)
    bam_to_pat(config)
    if not os.path.exists(config.get('blocks_file', '')):
        beta_to_blocks(config)
    else:
        os.system(f"cp {config['blocks_file']} {os.path.join(config['working_directory'], 'wgbstools/subsample.blocks')}")
    beta_to_table(config)
    coverage_to_illumina(config)
    illumina_to_methatlas(config)