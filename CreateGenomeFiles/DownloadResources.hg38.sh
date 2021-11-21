#set tool paths to your system's
SAMTOOLS=samtools-0.1.18
STAR=STAR252
BEDTOOLS=bedtools.2.27.1
BLAT=/home/alu/gabayo2/tools/blat/blat

mkdir ./processedTables
#---------------------------------------------- Genome----------------------------------------------
# genome fasta
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.chroms.tar.gz
tar -xzvf hg38.analysisSet.chroms.tar.gz
cat ./hg38.analysisSet.chroms/*.fa > ./processedTables/hg38_all.fa

#fasta index
${SAMTOOLS} faidx ./processedTables/hg38_all.fa
#STAR genome index
${STAR} --runMode genomeGenerate --genomeDir ./processedTables  --genomeFastaFiles ./processedTables/hg38_all.fa 

#---------------------------------------------- SNPs----------------------------------------------
# dbSNP tables (exchange for whatever tables you need)
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/snp150Common.txt.gz

#single genomic SNPs
zcat ./snp150Common.txt.gz | awk 'BEGIN{FS=OFS="\t"}{if (($4-$3==1) && ($11=="genomic") && ($12=="single")) print $2,$3,$4,$5,"0","+"}' > ./processedTables/dbSNP150_hg38_single_genomic.bed

#indels
zcat snp150Common.txt.gz | awk 'BEGIN{FS=OFS="\t"}{if (($11=="genomic") && (($12=="deletion") || ($12=="in-del") || ($12=="insertion"))) print $2,$3,$4,$5,"0","+"}' > ./processedTables/dbSNP150_hg38_deletion_in-del_insertion_genomic.bed
${BEDTOOLS} slop -i ./processedTables/dbSNP150_hg38_deletion_in-del_insertion_genomic.bed -b 50 -g ./processedTables/chrNameLength.txt > ./processedTables/dbSNP150_hg38_deletion_in-del_insertion_genomic_50nt.bed

# ORSHAY - PLEASE explain a little in this table and the matching 146 one
cat ./processedTables/dbSNP150_hg38_single_genomic.bed | awk -F '\t' '{printf("%s",$1); for (i=2;i<=NF;i++){if ($i=="") $i="NA"; printf("\t%s",$i);} printf "\n"}' | sort --parallel=24 -k1,1 -k2,2n -k3,3n > dbSNP150_hg38_single_genomic_sorted.bed 
cat dbSNP150_hg38_single_genomic_sorted.bed | cut -f 8 | tr "ACGTacgt" "TGCAtgca" | paste dbSNP150_hg38_single_genomic_sorted.bed /dev/stdin > dbSNP150_hg38_single_genomic_sorted_refByStrand.bed
cat dbSNP150_hg38_single_genomic_sorted_refByStrand.bed | grep -v -P "^chr.+_.+_alt\t" > dbSNP150_hg38_single_genomic_sorted_refByStrand_noAltChr.bed
cat | awk -F"\t" '$8 ~ /[ACGTacgt]/' dbSNP150_hg38_single_genomic_sorted_refByStrand_noAltChr.bed > dbSNP150_hg38_single_genomic_sorted_refByStrand_noAltChr_onlyACGT.bed
cat dbSNP150_hg38_single_genomic_sorted_refByStrand_noAltChr_onlyACGT.bed | grep -E '1000GENOMES|TOPMED' > dbSNP150_hg38_single_genomic_sorted_refByStrand_noAltChr_onlyACGT_TGP_TOPMED.bed
cat dbSNP150_hg38_single_genomic_sorted_refByStrand_noAltChr_onlyACGT_TGP_TOPMED.bed | cut -f 22 | tr "ACGT" "TGCA" | paste dbSNP150_hg38_single_genomic_sorted_refByStrand_noAltChr_onlyACGT_TGP_TOPMED.bed /dev/stdin > dbSNP150_hg38_single_genomic_sorted_refByStrand_noAltChr_onlyACGT_TGP_TOPMED_revAlleles.bed
cat dbSNP150_hg38_single_genomic_sorted_refByStrand_noAltChr_onlyACGT_TGP_TOPMED_revAlleles.bed | cut -f 6 | tr "+-" "-+" | paste dbSNP150_hg38_single_genomic_sorted_refByStrand_noAltChr_onlyACGT_TGP_TOPMED_revAlleles.bed /dev/stdin > dbSNP150_hg38_single_genomic_sorted_refByStrand_noAltChr_onlyACGT_TGP_TOPMED_revAlleles_revStrand.bed
cat dbSNP150_hg38_single_genomic_sorted_refByStrand_noAltChr_onlyACGT_TGP_TOPMED_revAlleles_revStrand.bed | awk -F'\t' '{if ($21>0) {strand=$6; if ($18 ~ /ObservedMismatch/) strand=$28; ref=$8; if (strand=="-") ref=$26; alleles=$22; n_ann_snp=split($9,a,"/"); if (n_ann_snp==$21 && $18 ~ /InconsistentAlleles/ && alleles !~ /N/) alleles=$27; alleles=substr(alleles,1,length(alleles)-1); split(alleles,arr_alleles,","); split($24,freq,",");max_allele="";max_freq=0;for (i in arr_alleles) {if (freq[i]>max_freq) {max_freq=freq[i];max_allele=arr_alleles[i]}} if (max_allele!=ref) print $0}}' > ./processedTables/dbSNP150_hg38_single_genomic_sorted_refByStrand_noAltChr_onlyACGT_TGP_TOPMED_revAlleles_revStrand_wrongRefAllele.bed


wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/snp147Common.txt.gz
zcat snp147Common.txt.gz | awk 'BEGIN{FS=OFS="\t"}{if (($4-$3==1) && ($11=="genomic") && ($12=="single")) print $2,$3,$4,$5,"0","+"}' > ./processedTables/dbSNP147_hg38_single_genomic.bed

zcat snp147Common.txt.gz | awk 'BEGIN{FS=OFS="\t"}{if (($11=="genomic") && (($12=="deletion") || ($12=="in-del") || ($12=="insertion"))) print $2,$3,$4,$5,"0","+"}' > ./processedTables/dbSNP147_hg38_deletion_in-del_insertion_genomic.bed
${BEDTOOLS} slop -i ./processedTables/dbSNP147_hg38_deletion_in-del_insertion_genomic.bed -b 50 -g ./processedTables/chrNameLength.txt > ./processedTables/dbSNP147_hg38_deletion_in-del_insertion_genomic_50nt.bed


wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/snp146Common.txt.gz
zcat snp146Common.txt.gz | awk 'BEGIN{FS=OFS="\t"}{if (($4-$3==1) && ($11=="genomic") &&($12=="single")) print $2,$3,$4,$5,$6,$7}' > dbSNP146_hg38_single_genomic.bed
cat dbSNP146_hg38_single_genomic.bed | awk -F '\t' '{printf("%s",$1); for (i=2;i<=NF;i++){if ($i=="") $i="NA"; printf("\t%s",$i);} printf "\n"}' | sort --parallel=24 -k1,1 -k2,2n -k3,3n > dbSNP146_hg38_single_genomic_sorted.bed
cat dbSNP146_hg38_single_genomic_sorted.bed | cut -f 8 | tr "ACGTacgt" "TGCAtgca" | paste dbSNP146_hg38_single_genomic_sorted.bed /dev/stdin > dbSNP146_hg38_single_genomic_sorted_refByStrand.bed
cat dbSNP146_hg38_single_genomic_sorted_refByStrand.bed | grep -v -P '^chr.+_.+_alt\t' > dbSNP146_hg38_single_genomic_sorted_refByStrand_noAltChr.bed
cat dbSNP146_hg38_single_genomic_sorted_refByStrand_noAltChr.bed | awk -F'\t' '$8 ~ /[ACGTacgt]/' > dbSNP146_hg38_single_genomic_sorted_refByStrand_noAltChr_onlyACGT.bed
cat dbSNP146_hg38_single_genomic_sorted_refByStrand_noAltChr_onlyACGT.bed | grep -E '1000GENOMES|TOPMED' > dbSNP146_hg38_single_genomic_sorted_refByStrand_noAltChr_onlyACGT_TGP_TOPMED.bed
cat dbSNP146_hg38_single_genomic_sorted_refByStrand_noAltChr_onlyACGT_TGP_TOPMED.bed | cut -f 22 | tr "ACGT" "TGCA" | paste dbSNP146_hg38_single_genomic_sorted_refByStrand_noAltChr_onlyACGT_TGP_TOPMED.bed /dev/stdin > dbSNP146_hg38_single_genomic_sorted_refByStrand_noAltChr_onlyACGT_TGP_TOPMED_revAlleles.bed
cat dbSNP146_hg38_single_genomic_sorted_refByStrand_noAltChr_onlyACGT_TGP_TOPMED_revAlleles.bed | cut -f 6 | tr "+-" "-+" | paste dbSNP146_hg38_single_genomic_sorted_refByStrand_noAltChr_onlyACGT_TGP_TOPMED_revAlleles.bed /dev/stdin > dbSNP146_hg38_single_genomic_sorted_refByStrand_noAltChr_onlyACGT_TGP_TOPMED_revAlleles_revStrand.bed
cat dbSNP146_hg38_single_genomic_sorted_refByStrand_noAltChr_onlyACGT_TGP_TOPMED_revAlleles_revStrand.bed | awk -F'\t' '{if ($21>0) {strand=$6; if ($18 ~ /ObservedMismatch/) strand=$28; ref=$8; if (strand=="-") ref=$26; alleles=$22; n_ann_snp=split($9,a,"/"); if (n_ann_snp==$21 && $18 ~ /InconsistentAlleles/ && alleles !~ /N/) alleles=$27; alleles=substr(alleles,1,length(alleles)-1); split(alleles,arr_alleles,","); split($24,freq,",");max_allele="";max_freq=0;for (i in arr_alleles) {if (freq[i]>max_freq) {max_freq=freq[i];max_allele=arr_alleles[i]}} if (max_allele!=ref) print $0}}' > ./processedTables/dbSNP146_hg38_single_genomic_sorted_refByStrand_noAltChr_onlyACGT_TGP_TOPMED_revAlleles_revStrand_wrongRefAllele.bed


#---------------------------------------------- ORSHAY PLEASE EXPLAIN HERE AND UNDER WHAT IS THIS----------------------------------------------
wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh38p7/known_issues/b150_GMAF_nonRefAllele_TGP.txt.gz
wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh38p7/known_issues/b150_GMAF_nonRefAllele_TOPMED.txt.gz

cat b150_GMAF_nonRefAllele_T*.txt | cut -f 1 | sort | uniq | awk '{print "rs"$0}' | grep -w -f /dev/stdin <(zcat snp150.txt.gz | cut -f 2-7) >b150_GMAF_nonRefAllele_TGP_TOPMED.bed
${BEDTOOLS} intersect -a b150_GMAF_nonRefAllele_TGP_TOPMED.bed -b <(zcat snp150.txt.gz | cut -f 2-7,11,12) -wa -wb -sorted | awk '$4==$10 && $3-$2==1 && $13=="genomic" && $14=="single"' | cut -f 1-6 > ./processedTables/b150_GMAF_nonRefAllele_TGP_TOPMED_single_genomic.bed

${BEDTOOLS} intersect -a b150_GMAF_nonRefAllele_TGP_TOPMED.bed -b <(zcat snp150.txt.gz | cut -f 2-7,11,12) -wa -wb -sorted | awk '$4==$10 && ($13=="genomic") && (($14=="deletion") || ($14=="in-del") || ($14=="insertion"))' | cut -f 1-6 > ./processedTables/b150_GMAF_nonRefAllele_TGP_TOPMED_deletion_in-del_insertion_genomic.bed


# ----------------------------------------------------- Repeats --------------------------------------------------------------------------------
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz
zcat rmsk.txt.gz | awk 'BEGIN{FS=OFS="\t"}{if ($13=="Alu") print $6,$7,$8,$11,"0",$10}' | sort -k1,1 -k2,2n > ./processedTables/Alu_hg38_sorted.bed
#One may add here more tables e.g.
zcat rmsk.txt.gz | awk 'BEGIN{FS=OFS="\t"}{if ($13=="Alu") print $0}' | cut -f 6-8| sort -k1,1 -k2,2n > ./processedTables/RRNA_hg38_sorted.bed

# ------------------------------------------------------ RefSeq -------------------------------------------------------------------------------
#Retrieval: RefSeq_Curated_NM_CDS.bed input file was retrieved using UCSC table browser (http://genome.ucsc.edu/cgi-bin/hgTables). We used the RefSeq Curated table from the NCBI RefSeq track and filtered only records with name=”NM_*”. The filtration results were retrieved in BED format requiring only coding exons in the final output file.
${BEDTOOLS} getfasta -bed RefSeq_Curated_NM_CDS.bed -fi ./processedTables/hg38_all.fa -s -bedOut 2> /dev/null | awk 'BEGIN{FS=OFS="\t"}{split($4,a,"_");print $1,$2,$3,a[1]"_"a[2],$5,$6,$7}' | grep -i -v "NN" | perl sequence_splicer_gene_coords.pl /dev/stdin RefSeq_Curated_NM_CDS_spliced_noNs.fa
perl sequence_chopper_with_coords.pl RefSeq_Curated_NM_CDS_spliced_noNs.fa 76 19 RefSeq_Curated_NM_CDS_spliced_noNs_chopped_76_19.fa > sequence_chopper_with_coords.log
${BLAT} ./processedTables/hg38_all.fa ./RefSeq_Curated_NM_CDS_spliced_noNs_chopped_76_19.fa -stepSize=5 -repMatch=2253 -minScore=0 -minIdentity=0 -out=pslx -noHead RefSeq_Curated_NM_CDS_spliced_noNs_chopped_76_19.pslx
cat RefSeq_Curated_NM_CDS_spliced_noNs_chopped_76_19.pslx | awk '$1>=68 && $2>0 && $5==0 {aln_blocks_sum=0; split($19,aln_blocks,","); for (i=1;i<=$18;i++) aln_blocks_sum+=aln_blocks[i]; if ($1/aln_blocks_sum >= 0.96) print $0}' | awk 'BEGIN{FS=OFS="\t"; rev["A"]="T";rev["C"]="G";rev["G"]="C";rev["T"]="A"}{aln_length=$1; aln_strand=$9; n=split($10,query_name,"_"); query_chr=query_name[4]; query_strand=query_name[3]; if (query_strand=="-") query_strand_flag=1; else query_strand_flag=0; j=1; for (i=5; i<=n; i++) {split(query_name[i],query_intervals,"-"); query_intervals[1]+=1; for (w=query_intervals[1+(1*query_strand_flag)]; w!=query_intervals[2-(1*query_strand_flag)]+(1-(2*query_strand_flag));w=w+(1-(2*query_strand_flag))){query_coords[j]=w;j++}} n_aln=$18; split($19,aln_blocks,","); split($20,query_blocks,","); split($21,target_blocks,","); split($22,query_seqs,",");  split($23,target_seqs,","); if (aln_strand=="-") aln_strand_flag=1; else aln_strand_flag=0; for (cur_block=n_aln-(n_aln+1)*aln_strand_flag;cur_block!=0+(n_aln+1)*aln_strand_flag; cur_block=cur_block-1+2*aln_strand_flag){cur_aln_length=aln_blocks[cur_block]; for (pos=1; pos<=cur_aln_length;pos++) {query_nuc=toupper(substr(query_seqs[cur_block],pos,1)); target_nuc=toupper(substr(target_seqs[cur_block],pos,1)); if (query_nuc!=target_nuc) {if (aln_strand=="-"){query_nuc=rev[query_nuc];target_nuc=rev[target_nuc]} query_pos=pos+query_blocks[cur_block]; query_coord=query_coords[query_pos+(j-1-2*query_pos+1)*aln_strand_flag]; target_pos=pos+target_blocks[cur_block]-1; print query_chr,query_coord-1,query_coord, query_nuc target_nuc,"0", query_strand, $14 "_" target_pos "_" target_pos+1 "_" aln_strand}}}}' > ./processedTables/RefSeq_Curated_NM_CDS_spliced_noNs_chopped_76_19_mms.bed
