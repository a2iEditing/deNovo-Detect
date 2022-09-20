import os
import shutil
import logging
import argparse
import logomaker
import pandas as pd
from glob import glob
import multiprocessing
import subprocess as sp
import matplotlib.pyplot as plt
from multiprocessing import Pool
import pybedtools.bedtool as bed
import csv

newLine = '\n'
tab = '\t'

# User arguments
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                 description='Reads the tissue file and run script for every donor in file')
parser.add_argument('-a', '--alignment', dest='align_dir', action='store', required=True,
                    help='Path to Alignment files')
parser.add_argument('-i', '--input', dest='region_file', action='store', required=True,
                    help='Path to a list of regions sites for analyze | Required file type: bed')
parser.add_argument('-b', '--blat', dest='blat_file', action='store', required=True,
                    help='Path to a file of mismatches to filter (produced from BLAT)| Required file type: bed')
parser.add_argument('-o', '--output', dest='output_dir', action='store', required=True, help='Path to output directory')
parser.add_argument('-l', '--log', dest='log_dir', action='store', default='', help='Path to a log file')
parser.add_argument('-g', '--genome', dest='genome', action='store', required=True,
                    help='Path to the genome of samples')
parser.add_argument('-pre', '--prefix', dest='pre', action='store', default='SRR',
                    help='Prefix of alignment files | Default prefix: SRR')
parser.add_argument('-chrom_length', dest='chrom_length', action='store', required=True,
                    help='List of all chromosomes length for the genome in use')
parser.add_argument('-intersect_list', dest='intersect_list', action='store', required=True,
                    help='List of bed files to intersect and filter')
parser.add_argument('-r', '--reads', dest='min_reads', action='store', default=100,
                    help='Minimum number of reads for a site')
parser.add_argument('-e', '--editing', dest='min_percentage', action='store', type=float, default=0.02,
                    help='Minimum editing percentage for a site')
parser.add_argument('-processes', dest='processes', action='store', default=20, help='Maximum number of processes')
parser.add_argument('--start_from', dest='start_from', action='store', type=int, default=1,
                    help='Step to start with(to skip steps, steps 1-3)')
parser.add_argument('--end_at', dest='end_at', action='store', type=int, default=3,
                    help='Step to end with(to stop at a step, steps 1-3)')
parser.add_argument('--strand_insensitive_blat', dest='strand_insensitive_blat', action='store_true',
                    help='If set, will use use strand-insensitive BLAT based filtering')

arguments = parser.parse_args()

# Remove '/' from the end of path if given
arguments.align_dir = arguments.align_dir.rstrip('/')
arguments.output_dir = arguments.output_dir.rstrip('/')

# Maximum number of cpu for using
max_processes = min(multiprocessing.cpu_count(), arguments.processes)

# Create log file
log_suggest_name = f'{arguments.output_dir.split("/")[-1]}_{arguments.genome.split("/")[-1].split(".")[0]}.log'
if arguments.log_dir == '':
    log_file_path = f'{arguments.output_dir}/Logs/{log_suggest_name}'
else:
    log_file_path = f'{arguments.log_dir}/{log_suggest_name}' if arguments.log_dir.endwith('/') else arguments.log_dir

arguments.log_dir = log_file_path
os.makedirs(os.path.dirname(arguments.log_dir), exist_ok=True)
logging.basicConfig(filename=arguments.log_dir, format='[%(asctime)s] %(levelname)-8s %(message)s',
                    datefmt='%y-%m-%d %H:%M:%S', level=logging.INFO, filemode='w')

logging.info('''Start running Or-Shay's pipeline\n''')

# Programs version
install_dir = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
program_dict = {line.split(':')[0].strip(): line.split(':')[1].strip() for line in
                open(f"{install_dir}/Programs_versions/Programs_versions_dictionary.txt")}
program_dict["install_dir"] = install_dir

# Default paths and arguments
genome_dir = arguments.genome.split('/')[:-1]
analysis_dir = f'{arguments.output_dir}/Analysis'
analysis_all_samples = f'{analysis_dir}/all_samples'
outputDataDir = f'{arguments.output_dir}/AllDataPairs'
outputStatsDir = f'{arguments.output_dir}/Statistics'
known_dir = f'{arguments.output_dir}/Known'
known_files_dir = f'{known_dir}/known_files'
known_all_samples = f'{known_files_dir}/all_samples'
sites_list_dir = f'{known_dir}/sites_for_known'
final_original_dir = f'{known_dir}/finale_lists/original'
final_filtered_dir = f'{known_dir}/finale_lists/filtered'

# Create all directories
os.makedirs(analysis_all_samples, exist_ok=True)
os.makedirs(outputDataDir, exist_ok=True)
os.makedirs(outputStatsDir, exist_ok=True)
os.makedirs(known_all_samples, exist_ok=True)
os.makedirs(sites_list_dir, exist_ok=True)
os.makedirs(final_original_dir, exist_ok=True)
os.makedirs(final_filtered_dir, exist_ok=True)

# Lists
bam_files = [x for x in glob(f'{arguments.align_dir}/{arguments.pre}*/*.bam')]
fileName = [os.path.basename(x) for x in glob(f'{arguments.align_dir}/{arguments.pre}*')]
clusterList = '0 50 100 200 400 800 1600 3200 6400 12800 25600 51200 102400'.split(' ')
mismatch_type = 'AC AG AT CA CG CT GA GC GT TA TC TG'.split(' ')

logging.info('User arguments:\n' + '\n'.join(f'{tab * 8}{k}: {v}' for k, v in vars(arguments).items()) + newLine)
logging.info('Programs versions:\n' + newLine.join(
    [f'{tab * 8}{key}: {value}' for key, value in program_dict.items()]) + newLine)
logging.info(f'{len(fileName)} files found on the directory: {arguments.align_dir}\n' + newLine.join(
    [f'{tab * 8}{file}' for file in fileName]) + newLine)


def step_1():
    logging.info('Starting step 1')
    with Pool(max_processes) as pool:
        pool.map(start_analyze, fileName)

    logging.info('Finished step 1 successfully!\n')


def start_analyze(sample):
    cur_aln_dir = f'{arguments.align_dir}/{sample}'
    cur_ident_dir = f'{analysis_dir}/{sample}'
    os.makedirs(cur_ident_dir, exist_ok=True)

    # Read-level detection
    #   read_overlap_md_tag_fixer.pl:
    #       THIS SCRIPT GETS SAM FILE IN WHICH OVERLAPPING BASES IN ONE MATE WERE SOFT-CLIPPED AND ORIGINAL CIGAR WAS KEPT IN A SPECIFIC BAM TAG.
    #       THE SCRIPT LOOKS FOR THE MD TAG AND FIXES IT TO MATCH TO THE NEW SOFT-CLIPPED CIGAR.
    #       OUTPUTS SAM FILE TO STDOUT.
    #
    #   mms_by_read_pos_v7.pl:
    #       This script gets sam file and output bed (bed6) file with all mismatches in the read level. Several filters are available: base quality, mapping quality, avoiding reads ends, avoiding bps around splice sites and ignoring duplicated reads.
    #       The script will output the mismatches relative to the positive DNA strand unless strand-specific data is given. In this case the read which include the RNA sequence should be specified.
    #       The bed output is as follows (tab delimited): chromosome, start coordinate (0-based),end coordinate(1-based), bed name, bed score: base quality, strand: mismatch strand.
    #       The name field (4th field in the output) is designed (semi-column ";" delimited): read name, read orientation (F or R), read "end" (1 or 2), mismatch type, mismatch position along the read sequence, mismatch position along the alignment,alignment length, read flag.
    #       INDEXES are not reported.
    #       MD tag must be compatible with CIGAR string (e.g. number of matched/deleted/inserted bps etc). Discrepancy might occur after applying tools that change the CIGAR on the input BAM/SAM file.
    cmd = r"""{samtools} view -h {cur_aln_dir}/{pre}*.bam | {perl} {install_dir}/Scripts/Step_1/read_overlap_md_tag_fixer.pl /dev/stdin | {perl} {install_dir}/Scripts/Step_1/mms_by_read_pos_v7.pl /dev/stdin {cur_ident_dir}/{pre}_byRead.bed /dev/stdout {cur_ident_dir}/{pre}_excludedReads.txt {cur_ident_dir}/{pre}_Aligned.sortedByCoord.out.mdup.clipOverCigar.filtered.sam 0 "" 30 33 5 5 4 5 5 5 "ALL" 0 2 3""".format(
        install_dir=program_dict["install_dir"],
        samtools=f'{program_dict["samtools"]}',
        perl=f'{program_dict["perl"]}',
        cur_aln_dir=f'{cur_aln_dir}',
        pre=f'{sample}', cur_ident_dir=f'{cur_ident_dir}')
    sp.check_output(cmd, shell=True, stderr=sp.PIPE)

    # Sites extraction to BED format:
    #   5th field contains: 
    #      * of reads supporting the mismatch
    #      * of reads supporting the mismatch on the plus strand
    cmd = r"""awk 'BEGIN{{FS=OFS="\t"}}{{split($4,a,";"); mms=$1 "\t" $2 "\t" $3 "\t" a[4]; seen[mms]++; if (a[2]=="F") plus[mms]++}}END{{for (i in seen) {{if (plus[i]=="") plus[i]=0; print i,seen[i]";"plus[i],"+"}}}}' {path}_byRead.bed | sort -k1,1 -k2,2n -k4,4 > {path}_all.bed""".format(
        path=f'{cur_ident_dir}/{sample}')
    sp.check_output(cmd, shell=True, stderr=sp.PIPE)

    cmd = f'bzip2 -f {cur_ident_dir}/{sample}_byRead.bed &'
    sp.check_output(cmd, shell=True, stderr=sp.PIPE)

    # Intersection with databases to include / exclude sites
    cur_file = f"{cur_ident_dir}/{sample}_all.bed"
    cur_file_pre = cur_file.replace('.bed', '')  # Removing the bed suffix

    for line in open(f'{arguments.intersect_list}', 'r').readlines():
        bed_name, bed_path, flags = line.split()
        cmd = f'{program_dict["bedtools"]} intersect -a {cur_file_pre}.bed -b {bed_path} {flags} > {cur_file_pre}_{bed_name}.bed'
        sp.check_output(cmd, shell=True, stderr=sp.PIPE)
        os.remove(f'{cur_file_pre}.bed')
        cur_file_pre = f'{cur_file_pre}_{bed_name}'

    # Intersection and strand correction with the user region file
    # Awk throws overlapping exons from '-' and '+' strand together
    cmd = r"""{bedtools} intersect -a {cur_file_pre}.bed -b {region_file} -wa -wb | cut -f 1-5,12 | awk 'BEGIN{{FS=OFS="\t"}}{{mms=$1 "\t" $2 "\t" $3 "\t" $4 "\t" $5; strand[mms]=strand[mms] "," $6}}END{{for (i in strand){{split(substr(strand[i],2),strand_arr,","); strand_flag=1;for (j in strand_arr){{if (strand_arr[1]!=strand_arr[j]){{strand_flag=0;break}}}} if (strand_flag==1) print i,strand_arr[1]}}}}' | sort -k1,1 -k2,2n -k4,4 | awk 'BEGIN{{FS=OFS="\t"; rev["AC"]="TG"; rev["AG"]="TC";rev["AT"]="TA";rev["CA"]="GT";rev["CG"]="GC";rev["CT"]="GA";rev["GA"]="CT";rev["GC"]="CG";rev["GT"]="CA";rev["TA"]="AT";rev["TC"]="AG";rev["TG"]="AC"}}{{if ($6=="-") $4=rev[$4]; print $1,$2,$3,$4,$5,$6}}' >  {cur_file_pre}_strandByCDS.bed""".format(
        bedtools=f'{program_dict["bedtools"]}', cur_file_pre=f'{cur_file_pre}', region_file=f'{arguments.region_file}')
    sp.check_output(cmd, shell=True, stderr=sp.PIPE)

    # Calculating depth
    cmd = r"""{samtools} mpileup -x -l {cur_file_pre}_strandByCDS.bed -BQ0 -d10000000 --ff SECONDARY,UNMAP,DUP {path}_Aligned.sortedByCoord.out.mdup.clipOverCigar.filtered.sam > {cur_file_pre}_strandByCDS.pu""".format(
        samtools=f'{program_dict["samtools"]}', path=f'{cur_ident_dir}/{sample}', cur_file_pre=f'{cur_file_pre}')
    sp.check_output(cmd, shell=True, stderr=sp.PIPE)

    # mpileup_analyser.pl:
    #   this script gets mpileup output file and returns a bed-like file with counting of reads with each nucleotide / with indels =
    cmd = r"""{perl} {install_dir}/Scripts/Step_1/mpileup_analyser.pl {cur_file_pre}_strandByCDS.pu /dev/stdout 30 33 | {bedtools} intersect -a {cur_file_pre}_strandByCDS.bed -b stdin -wa -wb | awk 'BEGIN{{FS=OFS="\t"}}{{print $1,$2,$3,$4,$5";"$29";"$17,$6}}' > {cur_file_pre}_sByCDS_wDPStrand.bed""".format(
        install_dir=program_dict["install_dir"], perl=f'{program_dict["perl"]}', bedtools=f'{program_dict["bedtools"]}',
        cur_file_pre=f'{cur_file_pre}')
    sp.check_output(cmd, shell=True, stderr=sp.PIPE)

    # Remove unnecessary files
    os.remove(f"{cur_ident_dir}/{sample}_Aligned.sortedByCoord.out.mdup.clipOverCigar.filtered.sam")
    os.remove(f"{cur_file_pre}_strandByCDS.bed")
    os.remove(f"{cur_file_pre}_strandByCDS.pu")

    logging.info(f'{tab}Sample {sample} finished step 1 successfully!')


def step_2():
    logging.info('Starting step 2')
    if not (os.path.isfile(f'{analysis_all_samples}/allSamples_withDepthPerStrand.bed.gz') or os.path.isfile(f'{analysis_all_samples}/allSamples_withDepthPerStrand.bed')):
        # Merging sites from all samples
        cmd = r"""{GAWK} 'BEGIN{{FS=OFS="\t";SUBSEP="\t"}}{{sites[$1,$2,$3,$4][ARGIND]=$5}}END{{for (i in sites) {{printf("%s",i); for (j=1; j<ARGC;j++) {{string="NA"; if (sites[i][j]!="") string=sites[i][j]; printf("\t%s",string)}} printf("\n");delete sites[i]}}}}' {analysis_dir}/{pre}*/*_all_*Strand.bed > {all_samples}/allSamples_withDepthPerStrand.bed""".format(
            GAWK=f'{program_dict["GAWK"]}', all_samples=f'{analysis_all_samples}', pre=f'{arguments.pre}',
            analysis_dir=f'{analysis_dir}')
    else:
        cmd = r"""gunzip {all_samples}/allSamples_withDepthPerStrand.bed.gz""".format(
            all_samples=f'{analysis_all_samples}')
    
    if not os.path.isfile(f'{analysis_all_samples}/allSamples_withDepthPerStrand.bed'):
        sp.check_output(cmd, shell=True, stderr=sp.PIPE)
    
    logging.info('Started statistical analysis and scoring')
   
    # Statistical analysis and scoring
    cmd = r"""{perl} {install_dir}/Scripts/Step_2/binom9.pl {all_samples}/allSamples_withDepthPerStrand.bed $(ls -1 {analysis_dir}/{pre}*/*_all*Strand.bed* | wc -l) | sort -k5gr > {all_samples}/allSamples_withDepthPerStrand_scoreAndPvals.txt""".format(
        install_dir=program_dict["install_dir"], perl=f'{program_dict["perl"]}', all_samples=f'{analysis_all_samples}',
        pre=f'{arguments.pre}', analysis_dir=f'{analysis_dir}')
    sp.check_output(cmd, shell=True, stderr=sp.PIPE)

    cmd = r"""gzip -f {all_samples}/allSamples_withDepthPerStrand.bed &""".format(all_samples=f'{analysis_all_samples}')
    sp.check_output(cmd, shell=True, stderr=sp.PIPE)

    # Filtering by scores
    logging.info('Filtering by score')

    cmd = r"""awk -F'\t' '$17 >= 0.25 && $5 >= 1' {all_samples}/allSamples_withDepthPerStrand_scoreAndPvals.txt > {all_samples}/allSamples_withDepthPerStrand_scoreAndPvals_filtered.txt""".format(
        all_samples=f'{analysis_all_samples}')
    sp.check_output(cmd, shell=True, stderr=sp.PIPE)

    # Printing pre-clustering summary
    cmd = r"""cat {all_samples}/allSamples_withDepthPerStrand_scoreAndPvals_filtered.txt | {bedtools} intersect -a stdin -b {region_file} -wa -wb | awk 'BEGIN{{FS=OFS=SUBSEP="\t"}}{{strand[$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17]=strand[$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17] "," $23}}END{{for (i in strand){{split(substr(strand[i],2),strand_arr,","); strand_flag=1;for (j in strand_arr){{if (strand_arr[1]!=strand_arr[j]){{strand_flag=0;break}}}} if (strand_flag==1) print i,strand_arr[1]}}}}' | sort -k1,1 -k2,2n -k4,4 > {all_samples}/allSamples_withDepthPerStrand_scoreAndPvals_filtered_withStrand.txt""".format(
        all_samples=f'{analysis_all_samples}', bedtools=f'{program_dict["bedtools"]}',
        region_file=f'{arguments.region_file}')
    sp.check_output(cmd, shell=True, stderr=sp.PIPE)

    # Filtering all the unknown chromosome from the list
    logging.info('Clear unknown chromosomes')

    clean_mess_chr = pd.read_csv(
        f'{analysis_all_samples}/allSamples_withDepthPerStrand_scoreAndPvals_filtered_withStrand.txt',
        header=None, sep='\t', dtype=str)
    clean_mess_chr = clean_mess_chr[~clean_mess_chr[0].str.contains(r'_(?!$)')]
    clean_mess_chr.to_csv(
        f'{analysis_all_samples}/allSamples_withDepthPerStrand_scoreAndPvals_filtered_withStrand_withoutUnknownChromosome.txt',
        index=False, header=None, sep='\t')

    # Filtering BLAT mismtaches
    if arguments.strand_insensitive_blat:
        logging.info("Filter BLAT mismatches - strand insensitive")
        cmd = r"""awk 'BEGIN{{arr["AC"]="AC";arr["AG"]="AG";arr["AT"]="AT";arr["CA"]="CA";arr["CG"]="CG";arr["CT"]="CT";arr["GA"]="CT";arr["GC"]="CG";arr["GT"]="CA";arr["TA"]="AT";arr["TC"]="AG";arr["TG"]="AC";}}{{OFS="\t";print $1"#"arr[$4],$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18}}' {all_samples}/allSamples_withDepthPerStrand_scoreAndPvals_filtered_withStrand_withoutUnknownChromosome.txt > {all_samples}/allSamples_withDepthPerStrand_scoreAndPvals_filtered_withStrand_withoutUnknownChromosome.MMCoord.bed""".format(
            all_samples=f'{analysis_all_samples}')
        sp.check_output(cmd, shell=True, stderr=sp.PIPE)

        cmd = r"""awk 'BEGIN{{arr["AC"]="AC";arr["AG"]="AG";arr["AT"]="AT";arr["CA"]="CA";arr["CG"]="CG";arr["CT"]="CT";arr["GA"]="CT";arr["GC"]="CG";arr["GT"]="CA";arr["TA"]="AT";arr["TC"]="AG";arr["TG"]="AC"}}{{OFS="\t"; print $1"#"arr[$4],$2,$3,$6}}' {blat_file}|{bedtools} intersect -a {all_samples}/allSamples_withDepthPerStrand_scoreAndPvals_filtered_withStrand_withoutUnknownChromosome.MMCoord.bed -b stdin -v -wa|awk -F "[\t|#]" '{{OFS="\t"; print $1,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19}}' > {all_samples}/allSamples_withDepthPerStrand_scoreAndPvals_filtered_withStrand_withoutUnknownChromosome.NoBLAT.txt""".format(
            all_samples=f'{analysis_all_samples}', bedtools=f'{program_dict["bedtools"]}',
            blat_file=f'{arguments.blat_file}')
        sp.check_output(cmd, shell=True, stderr=sp.PIPE)
    else:
        logging.info("Filter BLAT mismatches - strand sensitive")
        cmd = r"""awk '{{OFS="\t";print $1"#"$4,$2,$3,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18}}' {all_samples}/allSamples_withDepthPerStrand_scoreAndPvals_filtered_withStrand_withoutUnknownChromosome.txt > {all_samples}/allSamples_withDepthPerStrand_scoreAndPvals_filtered_withStrand_withoutUnknownChromosome.MMCoord.bed""".format(
            all_samples=f'{analysis_all_samples}')
        sp.check_output(cmd, shell=True, stderr=sp.PIPE)

        cmd = r"""awk '{{OFS="\t"; print $1"#"$4,$2,$3,$6}}' {blat_file}|{bedtools} intersect -a {all_samples}/allSamples_withDepthPerStrand_scoreAndPvals_filtered_withStrand_withoutUnknownChromosome.MMCoord.bed -b stdin -v -wa|awk -F "[\t|#]" '{{OFS="\t"; print $1,$3,$4,$2,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18}}' > {all_samples}/allSamples_withDepthPerStrand_scoreAndPvals_filtered_withStrand_withoutUnknownChromosome.NoBLAT.txt""".format(
            all_samples=f'{analysis_all_samples}', bedtools=f'{program_dict["bedtools"]}',
            blat_file=f'{arguments.blat_file}')
        sp.check_output(cmd, shell=True, stderr=sp.PIPE)

    # use annovar to get mRNA positions
    cmd = r"""awk 'BEGIN{{FS="\t";OFS=" "}}{{print $1,$3,$3,substr($4,1,1),substr($4,2,2),$1"_"$3"_"substr($4,1,1)substr($4,2,2)}}' {all_samples}/allSamples_withDepthPerStrand_scoreAndPvals_filtered_withStrand_withoutUnknownChromosome.NoBLAT.txt >{all_samples}/allSamples_withDepthPerStrand_scoreAndPvals_filtered_withStrand_withoutUnknownChromosome.forAnnovar.txt""".format(
        all_samples=f'{analysis_all_samples}')
    sp.check_output(cmd, shell=True, stderr=sp.PIPE)

    cmd = r"""{perl} {annovar}/table_annovar.pl  -buildver hg38 {all_samples}/allSamples_withDepthPerStrand_scoreAndPvals_filtered_withStrand_withoutUnknownChromosome.forAnnovar.txt -protocol ncbiRefSeq {annovarDB} -operation gx -nastring . -csvout -polish  --otherinfo""".format(
        perl=f'{program_dict["perl"]}',
        annovar=f'{program_dict["annovar"]}',
        annovarDB=f'{program_dict["annovarDB"]}',
        all_samples=f'{analysis_all_samples}')
    sp.check_output(cmd, shell=True, stderr=sp.PIPE)

    annotated_coords = f"{analysis_all_samples}/allSamples_withDepthPerStrand_scoreAndPvals_filtered_withStrand_withoutUnknownChromosome.forAnnovar.txt.hg38_multianno.csv"
    refseq_recs = []
    with open(annotated_coords) as annovar_res:
        reader = csv.DictReader(annovar_res, delimiter=",", lineterminator="\n")
        for coord in reader:
            refseqs = coord["AAChange.ncbiRefSeq"]
            if refseqs == ".":
                orig_pos = coord["Otherinfo1"]
                ref_id = coord["Gene.ncbiRefSeq"]
                refseq_recs.append("\t".join([ref_id, coord["Start"], coord["End"], orig_pos, "0", "+"]))
                continue
            refseqs = refseqs.strip(",").split(",")
            for refseq in refseqs:
                try:
                    _, ref_id, _, pos, _ = refseq.split(":")
                    pos = pos[3:-1]
                    orig_pos = coord["Otherinfo1"]
                    refseq_recs.append("\t".join([ref_id, str(int(pos) - 1), pos, orig_pos, "0", "+"]))
                except ValueError:
                    pass

    refeseqed_bed = f"{analysis_all_samples}/allSamples_withDepthPerStrand_scoreAndPvals_filtered_withStrand_withoutUnknownChromosome.ByRefseq.bed"

    with open(refeseqed_bed, 'w') as refeseqed_bed_o:
        refeseqed_bed_o.write("\n".join(refseq_recs))

    logging.info(f'{tab}Starting create clusters list')
    with Pool(min(len(clusterList), max_processes)) as cluster_pool:
        cluster_pool.map(clustering_files, clusterList)

    clustering_table = pd.DataFrame()
    clustering_table['Type'] = ['AC', 'AG', 'AT', 'CA', 'CG', 'CT', 'GA', 'GC', 'GT', 'TA', 'TC', 'TG', '%AG',
                                '%(AG+TC)', '%AG(!TC)']

    # Merging a clustering table from temp files
    for cluster_size in clusterList:
        temp_file = ''.join(glob(f"{outputDataDir}/*_{cluster_size}.tmp"))
        temp = pd.read_table(temp_file, dtype='str')
        clustering_table[temp.columns[0]] = temp[temp.columns[0]]
        os.remove(temp_file)

    clustering_table.to_csv(
        f'{outputStatsDir}/allSamples_withDepthPerStrand_scoreAndPvals_filtered_PairsClustStats.txt', index=False,
        sep='\t')

    # Create a logo graph for each miss and cluster size
    logging.info(f'{tab}Starting create Logo graphs')
    with Pool(min(len(clusterList), max_processes)) as logo_pool:
        logo_pool.map(create_logo_by_cluster_size, clusterList)

    logging.info('Finished step 2 successfully!\n')


def clustering_files(cluster_size):
    # Create file with sites coordinates foreach cluster.
    cluster_filter(cluster_size)

    # temp files for building stat table with all clusters
    cmd = r"""echo {cluster_size} > {outputDataDir}/allSamples_withDepthPerStrand_scoreAndPvals_filtered_clustStats_{cluster_size}.tmp""".format(
        cluster_size=f'{cluster_size}', outputDataDir=f'{outputDataDir}')
    sp.check_output(cmd, shell=True, stderr=sp.PIPE)

    cmd = r"""cat {outputDataDir}/allSamples_withDepthPerStrand_scoreAndPvals_filtered_clust_{cluster_size}.txt | {GAWK} 'BEGIN{{FS=OFS="\t";mms["AC"]=0; mms["AG"]=0; mms["AT"]=0; mms["CA"]=0; mms["CG"]=0; mms["CT"]=0; mms["GA"]=0; mms["GC"]=0; mms["GT"]=0; mms["TA"]=0; mms["TC"]=0; mms["TG"]=0}}{{mms[$4]++}}END{{n=asorti(mms,mms_sort);for (i=1; i<=n; i++) {{mms_type=mms_sort[i];mms_count=mms[mms_sort[i]];total+=mms_count; if (mms_type=="AG") ag=mms_count; if (mms_type=="TC") tc=mms_count; print mms[mms_sort[i]]}}print ag/total; print (ag+tc)/total; print ag/(total-tc)}}' >> {outputDataDir}/allSamples_withDepthPerStrand_scoreAndPvals_filtered_clustStats_{cluster_size}.tmp""".format(
        outputDataDir=f'{outputDataDir}', cluster_size=f'{cluster_size}', GAWK=f'{program_dict["GAWK"]}')
    sp.check_output(cmd, shell=True, stderr=sp.PIPE)


def cluster_filter(cluster_size):
    cmd = r"""{bedtools} sort -i {all_samples}/allSamples_withDepthPerStrand_scoreAndPvals_filtered_withStrand_withoutUnknownChromosome.ByRefseq.bed| {bedtools}  cluster -i stdin  -d {cluster_size}|cut -f 4,7|awk -F"[\t|_]" '{{OFS= ","; print $1,$2,$3,$4}}' > {all_samples}/allSamples_clust_{cluster_size}.ByRefseq.csv""".format(
        bedtools=f'{program_dict["bedtools"]}',
        all_samples=f'{analysis_all_samples}',
        cluster_size=f'{cluster_size}')
    sp.check_output(cmd, shell=True, stderr=sp.PIPE)

    df = pd.read_csv(f'{analysis_all_samples}/allSamples_clust_{cluster_size}.ByRefseq.csv',
                     names=["chr", "pos", "mutation_type", "cluster_num"],
                     sep=',')

    result_df = df.copy()
    result_df['marked_for_deletion'] = False

    cluster_size = int(cluster_size)
    for curr_index, curr_data in df.iterrows():
        curr_mutation = curr_data["mutation_type"]
        curr_cluster = curr_data["cluster_num"]

        cluster_df = df[(df["cluster_num"] == curr_cluster)
                        & (df["mutation_type"] != curr_mutation)]

        if len(cluster_df):
            result_df.loc[curr_index, "marked_for_deletion"] = True
            result_df.loc[cluster_df.index, 'marked_for_deletion'] = True
    result_df = result_df.drop("cluster_num", 1)
    result_df = result_df.drop_duplicates()
    any_refseq = pd.DataFrame(result_df.groupby(["chr", "pos"])["marked_for_deletion"].any())
    any_refseq = any_refseq.rename(columns={"marked_for_deletion": "to_delete"})
    result_df = result_df.merge(any_refseq, left_on=["chr", "pos"], right_index=True).drop("marked_for_deletion",
                                                                                           1).drop_duplicates()
    all_sites = pd.read_csv(
        f'{analysis_all_samples}/allSamples_withDepthPerStrand_scoreAndPvals_filtered_withStrand_withoutUnknownChromosome.NoBLAT.txt',
        sep='\t', header=None,
        names=["chr", "pre_pos", "pos", "mutation_type", "info_1", "info_2", "info_3", "info_4", "info_5",
               "info_6", "info_7", "info_8", "info_9", "info_10", "info_11", "info_12", "info_13", "strand"])
    result_df = result_df.merge(all_sites, on=["chr", "pos", "mutation_type"], sort=True)
    result_df = result_df[~result_df["to_delete"]]
    result_df.drop(['to_delete'], axis=1, inplace=True)
    result_df = result_df[["chr", "pre_pos", "pos", "mutation_type", "info_1", "info_2", "info_3", "info_4", "info_5",
                           "info_6", "info_7", "info_8", "info_9", "info_10", "info_11", "info_12", "info_13",
                           "strand"]]
    result_df.to_csv(f'{outputDataDir}/allSamples_withDepthPerStrand_scoreAndPvals_filtered_clust_{cluster_size}.txt',
                     header=None, index=False, sep='\t')

    logging.info(f'{tab * 2}Cluster size {cluster_size} finished successfully!')


def create_logo_by_cluster_size(cluster_size):
    # Create a Logo folder
    logo_dir = f'{arguments.output_dir}/Logo'

    # Iterating threw all mismatches types
    for miss in [x for x in mismatch_type]:

        # Create a folder for temp files
        temp_file_dir = f'{logo_dir}/{cluster_size}/temp_files'
        os.makedirs(temp_file_dir, exist_ok=True)

        # Read the file of certain cluster size and filter only the current mismatch
        miss_table = pd.read_csv(
            f'{outputDataDir}/allSamples_withDepthPerStrand_scoreAndPvals_filtered_clust_{cluster_size}.txt',
            usecols=[0, 1, 2, 3, 4, 17], sep='\t', header=None)
        miss_table = miss_table[miss_table[3] == miss]

        # Check if there is sites within the current mismatch
        if len(miss_table) > 1:
            # Convert from pandas format to Bedtools format
            miss_table = bed.BedTool.from_dataframe(miss_table)

            # Operate 'slop' function
            miss_table = miss_table.slop(g=arguments.chrom_length, b=2)

            # # Operate 'getFasta' function
            miss_table.sequence(fi=arguments.genome, fo=f'{temp_file_dir}/{miss}_getFasta.bed', s=True, tab=True)

            # Read only the sequence column from fasta file and save it as txt file
            get_fasta = pd.read_csv(f'{temp_file_dir}/{miss}_getFasta.bed', sep='\t', header=None, usecols=[1])
            get_fasta.to_csv(f'{temp_file_dir}/{miss}_getFasta_seqCols.txt', index=False, header=False)

            # Create a wight matrix for logo graph
            cmd = fr"""{program_dict["perl"]} {program_dict["install_dir"]}/Scripts/Step_2/logo_wight_matrix_calculator_dna.pl {temp_file_dir}/{miss}_getFasta_seqCols.txt {temp_file_dir}/{miss}_perl_out.txt -2"""
            sp.check_output(cmd, shell=True, stderr=sp.PIPE)

            # Create logo graph
            logo_graph = pd.read_csv(f'{temp_file_dir}/{miss}_perl_out.txt', sep='\t', index_col=0)
            crp_logo = logomaker.Logo(logo_graph, shade_below=.5, fade_below=.5, font_name='Arial Rounded MT Bold')
            crp_logo.style_spines(visible=False)
            crp_logo.style_spines(spines=('left', 'bottom'), visible=True)
            crp_logo.style_xticks(rotation=90, fmt='%d', anchor=0)

            plt.savefig(f'{logo_dir}/{cluster_size}/clust_{cluster_size}_miss_{miss}.pdf', dpi=600)
            logging.info(
                f'{tab * 2}Logo graphs of mismatch: {miss} for cluster size :{cluster_size} created successfully!')

        # Remove temp files folder
        shutil.rmtree(f'{temp_file_dir}/')


def step_3():
    logging.info('Staring step 3')
    # Aggregate all tissues sites for known for all clusters size
    for cluster in clusterList:
        cmd_aggregate = r'''cat {outputDataDir}/allSamples_withDepthPerStrand_scoreAndPvals_filtered_clust_{cluster}.txt | awk 'BEGIN{{FS="[\t]"; OFS="\t"}}{{print $1,$2,$3,$4,0,$18}}' | sort -k1,1 -k2,2n -k3,3n -k4 -u > {sites_list_dir}/sites_for_known_clust_{cluster}.bed'''.format(
            outputDataDir=f'{outputDataDir}', sites_list_dir=f'{sites_list_dir}', cluster=f'{cluster}')
        sp.check_output(cmd_aggregate, shell=True, stderr=sp.PIPE)

    logging.info(f'{tab}Staring known')
    with Pool(max_processes) as step_3_pool:
        step_3_pool.map(mpileup_known, bam_files)

    # Joining all samples
    cmd_joining = r'''{GAWK} 'BEGIN{{FS=OFS="\t";SUBSEP="\t"}}{{sites[$1,$2,$3,$4,"0",$6][ARGIND]=$5}}END{{for (i in sites) {{printf("%s",i); for (j=1; j<ARGC;j++) {{string="NA"; if (sites[i][j]!="") string=sites[i][j]; printf("\t%s",string)}} printf("\n");delete sites[i]}}}}'  $(ls -1 {known_files_dir}/*/*_all_*.bed) > {known_all_samples}/all_samples.txt'''.format(
        GAWK=f'{program_dict["GAWK"]}', known_files_dir=f'{known_files_dir}', known_all_samples=f'{known_all_samples}')
    sp.check_output(cmd_joining, shell=True, stderr=sp.PIPE)

    # Create file with sums and files of stats
    cmd_sum = r'''for file in {known_all_samples}/*all_samples.txt; do pre=${{file%.txt}}; awk -v p="$pre" 'BEGIN{{FS=OFS="\t";split("AC\tAG\tAT\tCA\tCG\tCT\tGA\tGC\tGT\tTA\tTC\tTG",types,"\t"); for (i in types) {{sup_ind[types[i]]=0;cov_ind[types[i]]=0}}}}{{if ($7 !~ /^[ACGTN]/) {{tot=0;sup=0;for (i=7;i<=NF;i++){{split ($i,a,";"); tot+=a[1];sup+=a[2];cov_ind[$4]+=a[1];sup_ind[$4]+=a[2]}} $5=tot ";" sup}} else $5="NA";print $0 > p"_wSums.txt"}}END{{for (i=1;i<=12;i++) {{print types[i],sup_ind[types[i]]/cov_ind[types[i]] > p"_indexStats.txt"}}}}' "$file"; done'''.format(
        known_all_samples=f'{known_all_samples}')
    sp.check_output(cmd_sum, shell=True, stderr=sp.PIPE)

    # Summarize to a final sites list
    cmd_w_sums_summary = r'''{GAWK} 'BEGIN{{FS=OFS="\t";SUBSEP="\t"}}{{sites[$1,$2,$3,$4,"0",$6][ARGIND]=$5}}END{{for (i in sites) {{printf("%s",i); for (j=1; j<ARGC;j++) {{string="NaN"; if (sites[i][j]!="") string=sites[i][j]; printf("\t%s",string)}} printf("\n");delete sites[i]}}}}' {known_all_samples}/all_samples_wSums.txt > {final_original_dir}/allTissues_wSums_summary_clust_0.txt'''.format(
        GAWK=f'{program_dict["GAWK"]}', known_all_samples=f'{known_all_samples}',
        final_original_dir=f'{final_original_dir}')
    sp.check_output(cmd_w_sums_summary, shell=True, stderr=sp.PIPE)

    columns = ['chr', 'start', 'end', 'mismatch', '0', 'strand', 'coverage', 'total coverage', 'edited coverage']
    cluster_0_list = pd.read_csv(f'{final_original_dir}/allTissues_wSums_summary_clust_0.txt', sep='\t', header=None)
    for cluster in clusterList:
        cluster_list = pd.read_csv(f'{sites_list_dir}/sites_for_known_clust_{cluster}.bed', sep='\t', header=None)
        cluster_list = pd.merge(left=cluster_0_list, right=cluster_list, how='right')
        cluster_list.to_csv(f'{final_original_dir}/allTissues_wSums_summary_clust_{cluster}.txt', header=False,
                            sep='\t', index=False)

        # Filtering minimum reads and editing percentage
        cluster_list = pd.concat([cluster_list, cluster_list.iloc[:, -1].str.split(';', expand=True)], axis=1)
        cluster_list.columns = columns
        cluster_list[['total coverage', 'edited coverage']] = cluster_list[['total coverage', 'edited coverage']].apply(
            pd.to_numeric)
        cluster_list = cluster_list.iloc[
                       ((cluster_list['edited coverage'] / cluster_list[
                           'total coverage']).values >= arguments.min_percentage) & (
                               cluster_list['total coverage'].values >= arguments.min_reads), :-2]
        cluster_list.to_csv(
            f'{final_filtered_dir}/allTissues_wSums_summary_clust_{cluster}_filtered_{arguments.min_reads}reads_and_{arguments.min_percentage}editingPercentage.txt',
            header=None, sep='\t', index=False)

    logging.info('Finished step 3 successfully!\n')


def mpileup_known(bam_file_path):
    known_pre = bam_file_path.split('/')[-2]
    sample_dir = f'{known_files_dir}/{known_pre}'
    os.makedirs(f'{sample_dir}', exist_ok=True)
    cmd_known = r'''{samtools} mpileup -l {sites_for_known} -BQ0 -d10000000 --ff SECONDARY,UNMAP,DUP {bam_file} | {perl} {install_dir}/Scripts/Step_3/mpileup_analyser.pl /dev/stdin {sample_dir}/{known_pre}_wDepth.txt 30 33'''.format(
        samtools=f'{program_dict["samtools"]}', sites_for_known=f'{sites_list_dir}/sites_for_known_clust_0.bed',
        bam_file=f'{bam_file_path}', perl=f'{program_dict["perl"]}', known_pre=f'{known_pre}',
        sample_dir=f'{sample_dir}', install_dir=program_dict["install_dir"])
    sp.check_output(cmd_known, shell=True, stderr=sp.PIPE)

    cmd_known = r'''{bedtools} intersect -a {sites_for_known} -b {sample_dir}/{known_pre}_wDepth.txt -wa -wb -nonamecheck | awk 'BEGIN{{FS=OFS="\t"; nuc_plus["A"]=0;nuc_plus["C"]=1;nuc_plus["G"]=2;nuc_plus["T"]=3; nuc_minus["A"]=3; nuc_minus["C"]=2; nuc_minus["G"]=1;nuc_minus["T"]=0}}{{mut=substr($4,2,1);offset=24; if ($6=="+") offset+=nuc_plus[mut]; else offset+=nuc_minus[mut];print $1,$2,$3,$4, $29 ";" $offset ,$6}}' > {sample_dir}/{known_pre}_covered_wDepth.bed'''.format(
        bedtools=f'{program_dict["bedtools"]}', sites_for_known=f'{sites_list_dir}/sites_for_known_clust_0.bed',
        known_pre=f'{known_pre}', sample_dir=f'{sample_dir}')
    sp.check_output(cmd_known, shell=True, stderr=sp.PIPE)

    cmd_known = r'''{bedtools} intersect -a {sites_for_known} -b {sample_dir}/{known_pre}_covered_wDepth.bed -v -nonamecheck | awk 'BEGIN{{FS=OFS="\t"}}{{print $1,$2,$3,$4, "0;0", $6}}' > {sample_dir}/{known_pre}_notCovered_wDepth.bed'''.format(
        bedtools=f'{program_dict["bedtools"]}', sites_for_known=f'{sites_list_dir}/sites_for_known_clust_0.bed',
        known_pre=f'{known_pre}', sample_dir=f'{sample_dir}')
    sp.check_output(cmd_known, shell=True, stderr=sp.PIPE)

    cmd_known = r'''cat {sample_dir}/{known_pre}_covered_wDepth.bed {sample_dir}/{known_pre}_notCovered_wDepth.bed | sort -k1,1 -k2,2n -k4 > {sample_dir}/{known_pre}_all_wDepth.bed'''.format(
        known_pre=f'{known_pre}', sample_dir=f'{sample_dir}')
    sp.check_output(cmd_known, shell=True, stderr=sp.PIPE)

    logging.info(f'{tab * 2}Sample {bam_file_path.split("/")[-1].split("_")[0]} Finished known successfully!')


if arguments.start_from == 1 and arguments.end_at >= 1:
    # Step 1
    step_1()

if arguments.start_from <= 2 and arguments.end_at >= 2:
    # Step 2
    step_2()

# Step 3
if arguments.end_at >= 3:
    step_3()

logging.info('''Finished De-Novo Pipeline successfully!''')
