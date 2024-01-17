#!/bin/bash
###################################################################################################
#                                                                                                 #
# ECPP 1.1.0-beta                                                                                 #
# Exon capture pipline for phylogentics                                                           #
# Adnan Moussalli (Museum Victoria, Australia)                                                    #
#                                                                                                 #
###################################################################################################
shopt -s expand_aliases
source /etc/profile

TEMP=`getopt -o m:p:s:a:r:t:e:f:x: --long argm:,argp:,args:,arga:,argr:,argt:,arge:,argf:,argx: -n 'ECCP_1.2.0.sh' -- "$@"`
eval set -- "$TEMP"
# extract options and their arguments into variables.
while true ; do
    case "$1" in
        -m|--argm)
            case "$2" in
                "") ARG_M='all' ; shift 2 ;;
                *) ARG_M=$2 ; shift 2 ;;
            esac ;;
        -p|--argp)
            case "$2" in
                "") shift 2 ;;
                *) ARG_P=$2 ; shift 2 ;;
            esac ;;
        -s|--args)
            case "$2" in
                "") shift 2 ;;
                *) ARG_S=$2 ; shift 2 ;;
            esac ;;
        -a|--arga)
            case "$2" in
                "") shift 2 ;;
                *) ARG_A=$2 ; shift 2 ;;
            esac ;;
        -r|--argr)
            case "$2" in
                "") shift 2 ;;
                *) ARG_R=$2 ; shift 2 ;;
            esac ;;
        -t|--argt)
            case "$2" in
                "") ARG_T='1' ; shift 2 ;;
                *) ARG_T=$2 ; shift 2 ;;
            esac ;;
        -e|--arge)
            case "$2" in
                 "") shift 2 ;;
                *) ARG_E=$2 ; shift 2 ;;
            esac ;;
        -f|--argf)
            case "$2" in
                 "") shift 2 ;;
                *) ARG_F=$2 ; shift 2 ;;
            esac ;;
        -x|--argx)
            case "$2" in
                "") shift 2 ;;
                *) ARG_X=$2 ; shift 2 ;;
            esac ;;
        --) shift ; break ;;
        *) echo "Internal error!" ; exit 1 ;;
    esac
done

###################################################################################################
#                                                                                                 #
# MODULE 1 (CLEAN) - DEDUPLICATE, SUBSET AND TRIM                                                 #
#                                                                                                 #
###################################################################################################
if [ "$ARG_M" = "all" ] || [ "$ARG_M" = "clean" ]; then
echo "MODULE $ARG_M"
echo -e "\tPE sample name list = $ARG_P"
echo -e "\tSE sample name list = $ARG_S"
echo -e "\tAmino acid reference = $ARG_A"
echo -e "\tPath to raw reads = $ARG_R"
echo -e "\tNumber of threads for multi-thread steps = $ARG_T\n"

###################### Prepare output subfolders ##################################################

mkdir -p $PWD/1_reads/LAST &&
mkdir -p $PWD/1_reads/TRIM &&
mkdir -p $PWD/1_reads/DR &&
mkdir -p $PWD/7_sum &&

# Read PE sample names into array
readarray -t P < $ARG_P
count=0
for i in "${P[@]}"
do
	(( count++ ))
	echo -e "Processing $count of ${#P[@]} - $i"  &&

	printf '%s\t' "$(date)" > $PWD/1_reads/$i.dupsum.PE.log && 
	printf '%s\t' $i >> $PWD/1_reads/$i.dupsum.PE.log && 
	printf '%s\t' "$(zgrep '\s1:' $ARG_R/$i\_R1.fastq | wc -l)" >> $PWD/1_reads/$i.dupsum.PE.log &&

	# subset reads to only those having a hit to reference
	echo -e "\tLASTAL subset paired reads." &&
	pigz -d $ARG_R/$i\_R1.fastq.gz $ARG_R/$i\_R2.fastq.gz &&

	cat $ARG_R/$i\_R1.fastq | paste - - - - | awk '{print ">"$1 ; print $3;}' > $PWD/1_reads/LAST/$i\_R1.fasta &&
	cat $ARG_R/$i\_R2.fastq | paste - - - - | awk '{print ">"$1 ; print $3;}' > $PWD/1_reads/LAST/$i\_R2.fasta &&
	lastdb -p $PWD/1_reads/LAST/REFERENCE $ARG_A &&
	lastal \
		-pPAM30 \
		-j3 \
		-F15 \
		-fTab \
		-P$ARG_T \
		-D1e+10 \
		$PWD/1_reads/LAST/REFERENCE $ARG_A $PWD/1_reads/LAST/$i\_R1.fasta > $PWD/1_reads/LAST/$i\.R1.LASTAL.TAB.txt &&
	grep -v '#' $PWD/1_reads/LAST/$i\.R1.LASTAL.TAB.txt | awk '{print $7}' > $PWD/1_reads/LAST/$i\.R1.LASTAL.TAB.HITS.txt &&
	lastal \
		-pPAM30 \
		-j3 \
		-F15 \
		-fTab \
		-P$ARG_T \
		-D1e+10 \
		$PWD/1_reads/LAST/REFERENCE $ARG_A $PWD/1_reads/LAST/$i\_R2.fasta > $PWD/1_reads/LAST/$i\.R2.LASTAL.TAB.txt &&
	grep -v '#' $PWD/1_reads/LAST/$i\.R2.LASTAL.TAB.txt | awk '{print $7}' > $PWD/1_reads/LAST/$i\.R2.LASTAL.TAB.HITS.txt &&
	cat $PWD/1_reads/LAST/$i\.R1.LASTAL.TAB.HITS.txt $PWD/1_reads/LAST/$i\.R2.LASTAL.TAB.HITS.txt |\
	sort | uniq | sed -e s/@//g > $PWD/1_reads/LAST/$i\.HITS.PE.txt &&
	seqtk subseq $ARG_R/$i\_R1.fastq $PWD/1_reads/LAST/$i\.HITS.PE.txt > $PWD/1_reads/LAST/$i\_R1.SUB.fastq &&
	seqtk subseq $ARG_R/$i\_R2.fastq $PWD/1_reads/LAST/$i\.HITS.PE.txt > $PWD/1_reads/LAST/$i\_R2.SUB.fastq &&
	
	echo -e "\tTRIM paired reads."
	java -jar /usr/local/bin/trimmomatic.jar PE -phred33 -threads $ARG_T \
		$PWD/1_reads/LAST/$i\_R1.SUB.fastq \
		$PWD/1_reads/LAST/$i\_R2.SUB.fastq \
		$PWD/1_reads/TRIM/$i\_R1.TRIMP.fastq.gz \
		$PWD/1_reads/TRIM/$i\_R1.TRIMU.fastq.gz \
		$PWD/1_reads/TRIM/$i\_R2.TRIMP.fastq.gz \
		$PWD/1_reads/TRIM/$i\_R2.TRIMU.fastq.gz \
		ILLUMINACLIP:/usr/local/share/triMomatic/TruSeq3-PE-2.fa:1:35:15 LEADING:3 TRAILING:3 SLIDINGWINDOW:8:20 MINLEN:40 2>&1 >/dev/null |\
	awk '/Input\ Read\ Pairs:/ {print $4,$7,$8,$12,$13,$17,$18,$20,$21}' OFS="\t" > $PWD/1_reads/$i.subtrimsumPE.log &&
	rm -R $PWD/1_reads/LAST/$i\_R1.SUB.fastq $PWD/1_reads/LAST/$i\_R2.SUB.fastq
	
	pigz -d $PWD/1_reads/TRIM/$i\_R1.TRIMP.fastq.gz $PWD/1_reads/TRIM/$i\_R2.TRIMP.fastq.gz &&
	echo "$PWD/1_reads/TRIM/$i""_R1.TRIMP.fastq" >$PWD/1_reads/fastuniqfilelist.txt && 
	echo "$PWD/1_reads/TRIM/$i""_R2.TRIMP.fastq" >>$PWD/1_reads/fastuniqfilelist.txt &&
	fastuniq -i $PWD/1_reads/fastuniqfilelist.txt -t q \
		-o $PWD/1_reads/DR/$i\_R1.DR.fastq \
		-p $PWD/1_reads/DR/$i\_R2.DR.fastq -c 0 && 
	grep '\s1:' $PWD/1_reads/DR/$i\_R1.DR.fastq | wc -l >> $PWD/1_reads/$i.dupsum.PE.log  &&
	pigz $PWD/1_reads/TRIM/$i\_R1.TRIMP.fastq $PWD/1_reads/TRIM/$i\_R2.TRIMP.fastq &&
	pigz $PWD/1_reads/DR/$i\_R1.DR.fastq $PWD/1_reads/DR/$i\_R2.DR.fastq &&


#	dedupe.sh in1=$PWD/1_reads/TRIM/$i\_R1.TRIMP.fastq.gz in2=$PWD/1_reads/TRIM/$i\_R2.TRIMP.fastq.gz out=$PWD/1_reads/DR/$i\.deduped.fastq.gz pigz=t unpigz=t overwrite=true 1>>$PWD/1_reads/dedup.log 2>&1
#	reformat.sh in=$PWD/1_reads/DR/$i\.deduped.fastq.gz out1=$PWD/1_reads/DR/$i\_R1.DR.fastq out2=$PWD/1_reads/DR/$i\_R2.DR.fastq pigz=t unpigz=t overwrite=true 1>>$PWD/1_reads/reformat.log 2>&1
#	grep '\s1:' $PWD/1_reads/DR/$i\_R1.DR.fastq | wc -l >> $PWD/1_reads/$i.dupsum.PE.log 
#	pigz $PWD/1_reads/DR/$i\_R1.DR.fastq $PWD/1_reads/DR/$i\_R2.DR.fastq

	pigz $ARG_R/$i\_R1.fastq &&
	pigz $ARG_R/$i\_R2.fastq &&
	paste $PWD/1_reads/$i.dupsum.PE.log $PWD/1_reads/$i.subtrimsumPE.log >> $PWD/1_reads/$i.dupsubtrimsumPE.log; rm $PWD/1_reads/$i.subtrimsumPE.log
done

# Read SE sample names into array
if [ "$ARG_S" ]; then
	readarray -t I < $ARG_S 
	for i in "${I[@]}"
	do
		echo "DEDUPLICATE single reads $i"
		printf '%s\t' "$(date)" > $PWD/1_reads/$i.dupsum.SE.log && 
		printf '%s\t' $i >> $PWD/1_reads/$i.dupsum.SE.log && 
		printf '%s\t' "$(zgrep '\s1:' $ARG_R/$i\_RS.fastq.gz | wc -l)" >> $PWD/1_reads/$i.dupsum.SE.log &&

		dedupe.sh in=$ARG_R/$i\_RS.fastq.gz out=$PWD/1_reads/DR/$i\_RS.DR.fastq pigz=t unpigz=t overwrite=true showspeed=t
		grep '\s1:' $PWD/1_reads/DR/$i\_RS.DR.fastq | wc -l >>$PWD/1_reads/$i.dupsum.SE.log &&

		echo "LASTAL subset single reads $i"	
		cat $PWD/1_reads/DR/$i\_RS.DR.fastq | paste - - - - | awk '{print ">"$1 ; print $3;}' > $PWD/1_reads/DR/$i\_RS.DR.fasta &&
		lastdb -p $PWD/1_reads/LAST/REFERENCE $ARG_A
		lastal \
			-pPAM30 \
			-j3 \
			-F15 \
			-fTab \
			-P$ARG_T \
			-D1e+3 \
			$PWD/1_reads/LAST/REFERENCE $ARG_A $PWD/1_reads/DR/$i\_RS.DR.fasta > $PWD/1_reads/LAST/$i\.RS.LASTAL.TAB.txt &&
		grep -v '#' $PWD/1_reads/LAST/$i\.RS.LASTAL.TAB.txt | awk '{print $7}' > $PWD/1_reads/LAST/$i\.RS.LASTAL.TAB.HITS.txt &&
		cat $PWD/1_reads/LAST/$i\.RS.LASTAL.TAB.HITS.txt |\
		sort | uniq | sed -e s/@//g > $PWD/1_reads/DR/$i\.HITS.SE.txt &&
		seqtk subseq $PWD/1_reads/DR/$i\_RS.DR.fastq $PWD/1_reads/DR/$i\.HITS.SE.txt > $PWD/1_reads/DRSUB/$i\_RS.DRSUB.fastq &&
	
		rm -R $PWD/1_reads/DR/$i\_RS.DR.fastq $PWD/1_reads/DR/$i\_RS.DR.fasta &&
		echo "TRIM single reads $i"
		java -jar /usr/local/bin/triMomatic.jar SE -phred33 -threads $ARG_T \
			$PWD/1_reads/DRSUB/$i\_RS.DRSUB.fastq \
			$PWD/1_reads/DRTRIM/$i\_RS.DRTRIMSE.fastq.gz \
			ILLUMINACLIP:/usr/local/share/triMomatic/TruSeq3-PE-2.fa:1:35:15 LEADING:3 TRAILING:3 SLIDINGWINDOW:8:20 MINLEN:40 2>&1 >/dev/null | \
		awk '/Input\ Read/ {print $3,$5,$6,$8,$9}' OFS="\t" > $PWD/1_reads/$i.subtrimsumSE.log &&
		rm -R $PWD/1_reads/DRSUB/$i\_RS.DRSUB.fastq &&
		paste $PWD/1_reads/$i.dupsum.SE.log $PWD/1_reads/$i.subtrimsumSE.log > $PWD/1_reads/$i.dupsubtrimsumSE.log
	done
fi
# Collate statistics
echo -e "Date\t Sample_Name\t RawReads\t DedupReads\t LastSubsetReads\t TrimPair\t TrimPair%\t TrimR1\t TrimR1%\t TrimR2\t TrimR2%\t TrimFail\t TrimFail%" > $PWD/7_sum/dupsubtrimsumPE.log &&
cat $PWD/1_reads/*.dupsubtrimsumPE.log | tr -d '()' >> $PWD/7_sum/dupsubtrimsumPE.log &&

if [ "$ARG_S" ]; then cat $PWD/1_reads/*.dupsubtrimsumSE.log > $PWD/7_sum/dupsubtrimsumSE.log; fi 

# Clean up
#rm $PWD/1_reads/*.dupsubtrimsumPE.log $PWD/1_reads/*.dupsum.PE.log 
#rm -rf $PWD/1_reads/LAST $PWD/1_reads/TRIM
fi	

###################################################################################################
#                                                                                                 #
# EMAP MODULE 2 (ASSEMBLE) - DE NOVO ASSEMBLY                                                     #
#                                                                                                 #
###################################################################################################
if [ "$ARG_M" = "all" ] || [ "$ARG_M" = "trinity" ] || [ "$ARG_M" = "spade" ]; then
echo "MODULE $ARG_M"
echo -e "\tPE sample name list = $ARG_P"
echo -e "\tSE sample name list = $ARG_S"
echo -e "\tAmino acid reference = $ARG_A"
echo -e "\tNumber of threads for multi-thread steps = $ARG_T\n"

###################### Prepare output subfolders ##################################################
mkdir -p $PWD/2_assemble &&
mkdir -p $PWD/3_blast &&
mkdir -p $PWD/3_blast/out &&
mkdir -p $PWD/4_map &&

###################### Read sample names into array ###############################################
readarray -t P < $ARG_P 
count=0
for i in "${P[@]}"
do
	# Spades Assembly #
	if [ "$ARG_M" = "spade" ]; then
	spades.py -k 21,33,55,77 --only-assembler --careful \
	-1 $PWD/1_reads/DRTRIM/$i\_R1.DRTRIMP.fastq.gz  \
	-2 $Rt.DRTRIMP.fastq \
	--s1 $Ft.DRTRIMU.fastq \
	--s2 $Rt.DRTRIMU.fastq \
	-o $PWD/assemble/$i -t $ARG_T

	# Make oneline fasta from spades assembly 
	# cat $PWD/assemble/$SN/contigs.fasta | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);} END {printf("\n");}' | awk 'NF' >$PWD/assemble/$i.spades.fasta && 
	# makesomethingNotInterleaved.pl $PWD/assemble/$SN/contigs.fasta > $PWD/assemble/$i.spades.fasta &&

	else

	# TRINITY Assembly #
	(( count++ ))
	echo -e "Assembling using TRINITY $count of ${#S[@]} - $i"

	# First concatenate paired and unpaired reads#
	if [ -f "$PWD/1_reads/DRTRIM/$i\_RS.DRTRIMSE.fastq.gz" ]; then
		cat 	$PWD/1_reads/DRTRIM/$i\_R1.DRTRIMP.fastq.gz \
			$PWD/1_reads/DRTRIM/$i\_R1.DRTRIMU.fastq.gz \
			$PWD/1_reads/DRTRIM/$i\_R2.DRTRIMU.fastq.gz \
			$PWD/1_reads/DRTRIM/$i\_RS.DRTRIMSE.fastq.gz > $PWD/1_reads/DRTRIM/$i\_R1.DRTRIMPU.fastq.gz
	else
		cat 	$PWD/1_reads/DR/$i\_R1.DR.fastq.gz \
			$PWD/1_reads/TRIM/$i\_R1.TRIMU.fastq.gz \
			$PWD/1_reads/TRIM/$i\_R2.TRIMU.fastq.gz > $PWD/1_reads/DR/$i\_R1.DRPU.fastq.gz
	fi
	Trinity --seqType fq --max_memory 40G --min_contig_length 100 --no_normalize_reads \
		--left $PWD/1_reads/DR/$i\_R1.DRPU.fastq.gz --right $PWD/1_reads/DR/$i\_R2.DR.fastq.gz  \
		--output $PWD/2_assemble/$i\.trinity --CPU $ARG_T --full_cleanup --min_kmer_cov 2 

#>>$PWD/2_assemble/Assembly.log 2>&1

	makesomethingNotInterleaved.pl $PWD/2_assemble/$i.trinity.Trinity.fasta > $PWD/2_assemble/$i.trinity.fasta; rm -f  $PWD/2_assemble/$i.trinity.Trinity.fasta
	fi

	echo "$i assemblies complete"
done
echo "ASSEMBLIES COMPLETE"
fi


###################################################################################################
#                                                                                                 #
# EMAP MODULE 3 (MAP) - MAPPING AND CONSENSUS CALLING                                             #
#                                                                                                 #
###################################################################################################
if [ "$ARG_M" = "all" ] || [ "$ARG_M" = "map" ]; then
echo "MODULE $ARG_M"
echo -e "\tPE sample name list = $ARG_P"
echo -e "\tSE sample name list = $ARG_S"

grep -e ">" $ARG_A |awk 'sub(/^>/, "")' > $PWD/7_sum/exonheaders.txt
sort -k1,1 $PWD/7_sum/exonheaders.txt | uniq > $PWD/7_sum/exonheader.sorted.uniq.txt

Nexons=$(wc -l < $PWD/7_sum/exonheader.sorted.uniq.txt) 
echo -e "\tTotal number of exons in the reference = $Nexons"
echo -e "\tNumber of threads for multi-thread steps = $ARG_T\n"

###################### Prepare output subfolders ##################################################
mkdir -p $PWD/4_map &&

###################### Read sample names into array ###############################################
readarray -t P < $ARG_P  
for i in "${P[@]}"
do
	## Make Blast Database and tBLASTn to aa reference (Can potentially replace this with LAST)
	makeblastdb -in $PWD/2_assemble/$i.trinity.fasta -title $i.trinity -out $PWD/3_blast/db/$i.trinity -dbtype nucl &&
	tblastn -query $ARG_A \
		-db $PWD/3_blast/db/$i.trinity -evalue 0.0001 \
		-seg no \
		-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen gaps ppos frames qseq" \
		-num_threads $ARG_T \
		> $PWD/3_blast/out/$i.tblastn &&

	### Sort and select top hit per reference exon.
	### The first sort orders the blast output by query exon name(ref exon) (-k1) then by bitscore (-k12) descending, then by evalue (-k11) ascending. 
	### The second sort picks the first (top) hit per query (ref exon) (i.e. the first best contig per reference exon).  
	sort -k1,1 -k12,12gr -k11,11g $PWD/3_blast/out/$i.tblastn | sort -u -k1,1 --merge >$PWD/3_blast/out/$i.tblastn.top &&

	### Sort by subject (contigs) and ensure each has only one hit.  Multiple hits are allowed if the exons come from the same gene.  If they do not, then pick only the best scoring.
	sort -k2,2 -k12,12gr -k11,11g $PWD/3_blast/out/$i.tblastn.top >$PWD/3_blast/out/$i.tblastn.top2 &&
	awk 'BEGIN { FS = OFS = "\t"} {split($1,a,"_"); $(NF+1)=a[1];$(NF+2)=$2; print}' $PWD/3_blast/out/$i.tblastn.top2 >$PWD/3_blast/out/$i.tblastn.top3 &&
	awk -F '\t' 'BEGIN {OFS=FS} {if ($2 == prev2 && $19==prev19) $2 = $2$1; else prev2 = $2; prev19=$19; print}' $PWD/3_blast/out/$i.tblastn.top3 | \
        sort -u -k2,2 --merge >$PWD/3_blast/out/$i.tblastn.top4 && 

	## Only keep contig if OHR is >50%, i.e. length of local hit is at least 50% of the reference length
	awk '{if ($4/$13>=0.01)print}' $PWD/3_blast/out/$i.tblastn.top4 >$PWD/3_blast/out/$i.tblastn.top5 &&

	## Only keep contig if mismatch is <30%, this is highly dependent on expected distance to reference
	awk '{if ($5/$4<=0.3)print}' $PWD/3_blast/out/$i.tblastn.top5 >$PWD/3_blast/out/$i.tblastn.top6 && 

	## Extract contigs using local alignment cordinates and creating a sample specific reference for mapping.  These will be in frame given tblastn was used. 
	pullexons_EC_AM2.py $PWD/2_assemble/$i.trinity.fasta $PWD/3_blast/out/$i.tblastn.top6 $PWD/4_map/$i.tblastn.top6.exons &&

	# Concatenate unpaired reads #
	if [ -f "$PWD/1_reads/DRTRIM/$i\_RS.DRTRIMSE.fastq.gz" ]; then
		cat 	$PWD/1_reads/DRTRIM/$i\_R1.DRTRIMU.fastq.gz \
			$PWD/1_reads/DRTRIM/$i\_R2.DRTRIMU.fastq.gz \
			$PWD/1_reads/DRTRIM/$i\_RS.DRTRIMSE.fastq.gz > $PWD/1_reads/DRTRIM/$i\_R1.DRTRIMUP.fastq.gz
	else
		cat 	$PWD/1_reads/TRIM/$i\_R1.TRIMU.fastq.gz \
			$PWD/1_reads/TRIM/$i\_R2.TRIMU.fastq.gz > $PWD/1_reads/TRIM/$i\_R1.TRIMUP.fastq.gz
	fi
echo "a"	
	bbwrap.sh threads=54 k=12 unpigz=t minratio=0.25 maxindel=100 local\
	ref=$PWD/4_map/$i.tblastn.top6.exons \
	in1=$PWD/1_reads/DR/$i\_R1.DR.fastq.gz,$PWD/1_reads/TRIM/$i\_R1.TRIMUP.fastq.gz \
	in2=$PWD/1_reads/DR/$i\_R2.DR.fastq.gz,null out=$PWD/4_map/$i.trinity_bbmap.sam  \
	append nodisk \
	statsfile=$PWD/4_map/$i.trinity_bbmap.stat covstats=$PWD/4_map/$i.trinity_bbmap.covstat &&
echo "b"
	echo "$i" > $PWD/4_map/$i.trinity_bbmap.stat.sum
	awk '/Reads\ Used:/ {print $3}' OFS="\t"  $PWD/4_map/$i.trinity_bbmap.stat >> $PWD/4_map/$i.trinity_bbmap.stat.sum
	awk '/mated\ pairs:/ {print $3}' OFS="\t"  $PWD/4_map/$i.trinity_bbmap.stat >> $PWD/4_map/$i.trinity_bbmap.stat.sum
	awk '/mapped:/ {print $2}' OFS="\t"  $PWD/4_map/$i.trinity_bbmap.stat >> $PWD/4_map/$i.trinity_bbmap.stat.sum

	join -t $'\t' -a 1 -1 1 -2 1 -e NULL -o 2.2 \
	<(sort -k1,1 $PWD/7_sum/exonheaders.txt) <(sort -k1,1 $PWD/4_map/$i.trinity_bbmap.covstat) \
	> $PWD/7_sum/$i.trinity_bbmap.covstat.avgfold &&
	echo -e "$i" | cat - $PWD/7_sum/$i.trinity_bbmap.covstat.avgfold > $PWD/4_map/$i.trinity_bbmap.covstat.avgfoldh; rm -f $PWD/7_sum/$i.trinity_bbmap.covstat.avgfold &&
	
	samtools view -S $PWD/4_map/$i.trinity_bbmap.sam -b -o $PWD/4_map/$i.trinity_bbmap.bam; rm -f $PWD/4_map/$i.trinity_bbmap.sam &&
	samtools sort -@ 54 $PWD/4_map/$i.trinity_bbmap.bam -o $PWD/4_map/$i.trinity_bbmap.sorted.bam; rm -f $PWD/4_map/$i.trinity_bbmap.bam &&
	samtools index $PWD/4_map/$i.trinity_bbmap.sorted.bam &&
	samtools mpileup -B -A -f $PWD/4_map/$i.tblastn.top6.exons $PWD/4_map/$i.trinity_bbmap.sorted.bam | \
	java -jar /usr/local/bin/varscan.jar mpileup2cns --min-coverage 4 -output-vcf 1 >$PWD/4_map/$i.trinity_bbmap_vscns.vcf && 
	pigz -f $PWD/4_map/$i.trinity_bbmap_vscns.vcf &&
	zcat $PWD/4_map/$i.trinity_bbmap_vscns.vcf.gz | vcf-to-tab > $PWD/4_map/$i.trinity_vscnss.tab &&
	samtools faidx $PWD/4_map/$i.tblastn.top6.exons &&
	fai2tab.pl $PWD/4_map/$i.tblastn.top6.exons.fai >$PWD/4_map/$i.tblastn.top6.exons.fai.tab &&
# the following awk works fine but seems to hang on very large (i.e. whole exome) data.  Not sure why.
	awk 'NR==FNR{a[$1,$2]=$3OFS$4;next}{$3=a[$1,$2];print}' OFS='\t' $PWD/4_map/$i.trinity_vscnss.tab $PWD/4_map/$i.tblastn.top6.exons.fai.tab >$PWD/4_map/$i.trinity_vscnsf.tab &&
#so used join with a bit of fiddle to get the same results - much faster in general anyway.
#	join -t $'\t' -a 1 -1 1 -2 1 -e NULL -o 1.2,1.3,2.4,2.5 \#
#	<(<$PWD/4_map/$i.tblastn.top6.exons.fai.tab awk '{print $1"-"$2"\t"$0}' | sort -k1b,1) \
#	<(< $PWD/4_map/$i.trinity_vscnss.tab awk '{print $1"-"$2"\t"$0}' | sort -k1b,1) | sort -k1,1 -k2,2n > $PWD/4_map/$i.trinity_vscnsf.tab &&

	echo -e "CHROM\tPOS\tREF\t$i" | cat - $PWD/4_map/$i.trinity_vscnsf.tab  >$PWD/4_map/$i.trinity_vscnsfh.tab &&
	awk 'BEGIN { FS = OFS = "\t" } { for(j=1; j<=NF; j++) if($j ~ /^ *$/) $j = "n\tn/n" }; 1' $PWD/4_map/$i.trinity_vscnsfh.tab >$PWD/4_map/$i.trinity_vscnsfhn.tab &&
	vcf_tab_to_fasta_alignmentAM2.pl $PWD/4_map/$i.trinity_vscnsfhn.tab >$PWD/4_map/$i.trinity_vscnsfhn.fasta &&
	echo $i mapped and consensus created
done
fi
if [ "$ARG_M" = "all" ] || [ "$ARG_M" = "map" ]; then
	echo -e "EXONID" | cat - $PWD/7_sum/exonheader.sorted.uniq.txt > $PWD/7_sum/exonheadersh.txt; #rm -f $PWD/7_sum/exonheadersorted.txt &&
	paste $PWD/4_map/*.stat.sum > $PWD/7_sum/statsumMatrix.txt &&
	paste $PWD/4_map/*.avgfoldh > $PWD/7_sum/avgfoldmatrix.txt &&
	paste $PWD/7_sum/exonheadersh.txt $PWD/7_sum/avgfoldmatrix.txt > $PWD/7_sum/CoverageMatrix.txt && rm $PWD/7_sum/avgfoldmatrix.txt &&
	echo "MAPPING COMPLETE"
fi

###################################################################################################
#                                                                                                 #
# EMAP MODULE 4 (ALIGN) - TRANSLATE, ALIGN WITH MAFFT, BACKTRANSLATE                              #
#                                                                                                 #
###################################################################################################
if [ "$ARG_M" = "all" ] || [ "$ARG_M" = "align" ]; then
echo "MODULE $ARG_M"
echo -e "\tPE sample name list = $ARG_P"
echo -e "\tSE sample name list = $ARG_S"

sort $PWD/7_sum/exonheaders.txt| uniq >$PWD/7_sum/exonheaders.unique.txt
Nexons=$(wc -l < $PWD/7_sum/exonheaders.unique.txt) 

echo -e "\tTotal number of exons in the reference = $Nexons"
echo -e "\tNumber of threads for multi-thread steps = $ARG_T\n"

###################### Prepare output subfolders ##################################################
Ntaxa=$(wc -l < "$ARG_P")

mkdir -p $PWD/5_align.T$Ntaxa &&

echo "MODULE $ARG_M"
echo "sample name list = $ARG_P"
echo "amino acid reference = $ARG_A"
echo "number of threads for multi-thread steps = $ARG_T"

esweep \
	-i $PWD/4_map \
	-o $PWD/5_align.T$Ntaxa \
	-t $PWD/$ARG_P \
	-e $PWD/7_sum/exonheaders.unique.txt \
	-f %s.tblastn.top6.exons.fai \
	-b %s.trinity_vscnsfhn.fasta \
	-c \
	-x 3\
	-v && 

#esweep_linux_v2 \
#	-i $PWD/4_map \
#	-o $PWD/5_align \
#	-t $PWD/$ARG_S \
#	-e $PWD/7_sum/exonheaders.txt \
#	-f %s.tblastn.top6.exons.fai \
#	-b %s.trinity_vscnsfhn.fasta \
#	-p Bombyx_mori.ASM15162v1.37.exons.N2.fai \
#	-s SpeciesList.transcriptome.47.txt \
#	-g %s.T47.G12498.03.50.15.1.0.fasta \
#	-c \
#	-x 3 \
#	-v && 




mv $PWD/5_align.T$Ntaxa/matrix_length.csv  $PWD/5_align.T$Ntaxa/matrix_lengthT.$Ntaxa.E$Nexons.csv

cat $PWD/7_sum/exonheaders.unique.txt | while read exons
do
	echo "Processing $PWD/5_align.T$Ntaxa/${exons}_aligned.fasta"
	perl -i -nle "s/-|~|n/N/g unless $. & 1; print " $PWD/5_align.T$Ntaxa/${exons}_aligned.fasta && 
	degapcodon2.pl $PWD/5_align.T$Ntaxa/${exons}_aligned.fasta > $PWD/5_align.T$Ntaxa/${exons}_alignedDG.fasta &&
	alignAM.py \
		-prot_outfile $PWD/5_align.T$Ntaxa/${exons}_mafft_aa.fasta \
		-nuc_outfile $PWD/5_align.T$Ntaxa/${exons}_mafft_nt.fasta \
		-aligner mafft -options " --quiet --localpair --maxiterate 1000 --reorder --anysymbol --thread -1 " \
		$PWD/5_align.T$Ntaxa/${exons}_alignedDG.fasta ; #rm $PWD/5_align.T$Ntaxa/${exons}_unalignedDG.fasta $PWD/5_align.T$Ntaxa/${exons}_unaligned.fasta &&

	perl -i -nle "s/n|N/-/g unless $. & 1; print " $PWD/5_align.T$Ntaxa/${exons}_mafft_nt.fasta && 
	selectSites.pl -x 3 $PWD/5_align.T$Ntaxa/${exons}_mafft_nt.fasta > $PWD/5_align.T$Ntaxa/${exons}_mafft_nt_trim.fasta; rm $PWD/5_align.T$Ntaxa/${exons}_mafft_nt.fasta $PWD/5_align.T$Ntaxa/${exons}_mafft_aa.fasta
	makesomethingNotInterleaved.pl $PWD/5_align.T$Ntaxa/${exons}_mafft_nt_trim.fasta > $PWD/5_align.T$Ntaxa/${exons}_mafft_nt_trimNI.fasta; rm $PWD/5_align.T$Ntaxa/${exons}_mafft_nt_trim.fasta
done

echo "FIRST PASS ALIGNMENTS COMPLETE"
echo -e "\tPlease now assess exons_length.csv and exons_hetero.csv in 7_sum folder"
echo -e "\tDecide on which subset of taxa and which subset of exons should proceed with to the final module" 
echo -e "\tUse -e and -fto bring in new lists of subset exons and taxa names"

fi

############################# MODULE5  SUBSET #####################################################
#                                                                                                 #
# New lists of exons and taxa names are to be provided here via the -x and -t switch              #
#                                                                                                 #
###################################################################################################
if [ "$ARG_M" = "all" ] || [ "$ARG_M" = "subset" ]; then

###################### Prepare output subfolders ##################################################
Ntaxa=$(wc -l < "$ARG_P")
subexon=$(wc -l < "$ARG_E")
subtaxa=$(wc -l < "$ARG_F") 

echo $Ntaxa
echo $subexon
echo $subtaxa

DATE=$(date +%Y%m%d) &&
mkdir -p $PWD/6_Final_Alignments.T$subtaxa.E$subexon &&

echo "MODULE $ARG_M"
echo "number of threads for multi-thread steps = $ARG_T"
echo "subset sample name list = $ARG_F"
echo "subset exon name list = $ARG_E"

esweep \
	-i $PWD/5_align.T$Ntaxa \
	-o $PWD/6_Final_Alignments.T$subtaxa.E$subexon \
	-n $PWD/$ARG_F \
	-k $PWD/$ARG_E \
	-d %s_mafft_nt_trimNI.fasta \
	-x \
	-v
#esweep_linux_v2 \
#	-i $PWD/4_map \
#	-o $PWD/5_align \
#	-t $PWD/$ARG_S \
#	-e $PWD/7_sum/exonheaders.txt \
#	-f %s.tblastn.top6.exons.fai \
#	-b %s.trinity_vscnsfhn.fasta \
#	-p Bombyx_mori.ASM15162v1.37.exons.N2.fai \
#	-s SpeciesList.transcriptome.47.txt \
#	-g %s.T47.G12498.03.50.15.1.0.fasta \
#	-c \
#	-x 3 \
#	-v && 

readarray ExonSubset < $ARG_E
for i in "${ExonSubset[@]}"
do
	EN=( $i )
	wc -L $PWD/5_align.T$Ntaxa/${EN}_mafft_nt_trimNI.fasta > $PWD/6_Final_Alignments.T$subtaxa.E$subexon/${EN}.len1.txt
						
####### remove duMy sequences arising from esweep before doing BMGE
        degapcodon2.pl $PWD/6_Final_Alignments.T$subtaxa.E$subexon/${EN}_unaligned.fasta > $PWD/6_Final_Alignments.T$subtaxa.E$subexon/${EN}_unaligneddg.fasta && rm $PWD/6_Final_Alignments.T$subtaxa.E$subexon/${EN}_unaligned.fasta
	perl -p -e "s{^}{TTT} if $. %2 ==0" $PWD/6_Final_Alignments.T$subtaxa.E$subexon/${EN}_unaligneddg.fasta \
		> $PWD/6_Final_Alignments.T$subtaxa.E$subexon/${EN}_unalignedt.fasta && rm $PWD/6_Final_Alignments.T$subtaxa.E$subexon/${EN}_unaligneddg.fasta
	selectSites.pl -x 3  $PWD/6_Final_Alignments.T$subtaxa.E$subexon/${EN}_unalignedt.fasta \
		>  $PWD/6_Final_Alignments.T$subtaxa.E$subexon/${EN}_unalignedtt.fasta && rm $PWD/6_Final_Alignments.T$subtaxa.E$subexon/${EN}_unalignedt.fasta
	makesomethingNotInterleaved.pl $PWD/6_Final_Alignments.T$subtaxa.E$subexon/${EN}_unalignedtt.fasta \
		> $PWD/6_Final_Alignments.T$subtaxa.E$subexon/${EN}_unalignedNI.fasta && rm $PWD/6_Final_Alignments.T$subtaxa.E$subexon/${EN}_unalignedtt.fasta 
#	perl -i -nle "s/-/N/g unless $. & 1; print " $PWD/6_Final_Alignments.T$subtaxa.E$subexon/${EN}_unalignedt.fasta &&

# Alignment cleanup using BMG with -g gap rate cut-off per column = 50%, and -h entropy score cut-off of 0.5 (less is more stringent) - this may not be neccessary (see manual), and -w window of 1bp
# The trimming is progressively more stringent as the BLOSUM index increases (e.g. BLOSUM95); reciprocally, the trimming is progressively more relaxed as the BLOSUM index is lower (e.g. BLOSUM30). 
# In practice, it is recommended to use BLOSUM95 with closely related sequences, and BLOSUM30 with distantly related sequences. PAM30 is a good all round - but better still would be to provide this as a switch.
	java -jar /usr/local/bin/BMGE100.jar \
		-i $PWD/6_Final_Alignments.T$subtaxa.E$subexon/${EN}_unalignedNI.fasta \
		-t CODON \
		-m BLOSUM62 \
		-h 0.5 \
		-w 1 \
		-g  0.5 \
		-oco123f $PWD/6_Final_Alignments.T$subtaxa.E$subexon/${EN}_aligned_mafft_nt.BMGE.fasta &&
	rm $PWD/6_Final_Alignments.T$subtaxa.E$subexon/${EN}_unalignedNI.fasta
	wc -L $PWD/6_Final_Alignments.T$subtaxa.E$subexon/${EN}_aligned_mafft_nt.BMGE.fasta > $PWD/6_Final_Alignments.T$subtaxa.E$subexon/${EN}.len2.txt
	selectSites.pl -s '4-' $PWD/6_Final_Alignments.T$subtaxa.E$subexon/${EN}_aligned_mafft_nt.BMGE.fasta \
		> $PWD/6_Final_Alignments.T$subtaxa.E$subexon/${EN}_aligned_mafft_nt.BMGEt.fasta && rm $PWD/6_Final_Alignments.T$subtaxa.E$subexon/${EN}_aligned_mafft_nt.BMGE.fasta
	makesomethingNotInterleaved.pl $PWD/6_Final_Alignments.T$subtaxa.E$subexon/${EN}_aligned_mafft_nt.BMGEt.fasta \
		> $PWD/6_Final_Alignments.T$subtaxa.E$subexon/${EN}_aligned_mafft_nt.BMGEtIN.fasta && rm $PWD/6_Final_Alignments.T$subtaxa.E$subexon/${EN}_aligned_mafft_nt.BMGEt.fasta
	degapcodon2.pl $PWD/6_Final_Alignments.T$subtaxa.E$subexon/${EN}_aligned_mafft_nt.BMGEtIN.fasta \
		> $PWD/6_Final_Alignments.T$subtaxa.E$subexon/${EN}_aligned_mafft_nt.BMGE.30.fasta && rm $PWD/6_Final_Alignments.T$subtaxa.E$subexon/${EN}_aligned_mafft_nt.BMGEtIN.fasta

done
######## REPEAT ESWEEP HERE SUCH THAT THE FINAL COMBINED.NEX PARTITION INFORMATION IS UPDATED.
esweep \
	-i $PWD/6_Final_Alignments.T$subtaxa.E$subexon/ \
	-o $PWD/6_Final_Alignments.T$subtaxa.E$subexon/ \
	-n $PWD/$ARG_F \
	-k $PWD/$ARG_E \
	-d %s_aligned_mafft_nt.BMGE.30.fasta \
	-x \
	-v
for i in "${ExonSubset[@]}"
do
	EN=( $i )
	mv $PWD/6_Final_Alignments.T$subtaxa.E$subexon/${EN}_unaligned.fasta $PWD/6_Final_Alignments.T$subtaxa.E$subexon/${EN}_mafft_BMGE_H07W10G05.fasta && rm $PWD/6_Final_Alignments.T$subtaxa.E$subexon/${EN}_aligned_mafft_nt.BMGE.30.fasta
done
cat $PWD/6_Final_Alignments.T$subtaxa.E$subexon/*len1.txt > $PWD/6_Final_Alignments.T$subtaxa.E$subexon/Alignmentlength1.txt && rm $PWD/6_Final_Alignments.T$subtaxa.E$subexon/*len1.txt &&
cat $PWD/6_Final_Alignments.T$subtaxa.E$subexon/*len2.txt > $PWD/6_Final_Alignments.T$subtaxa.E$subexon/Alignmentlength2.txt && rm $PWD/6_Final_Alignments.T$subtaxa.E$subexon/*len2.txt
echo "FINAL ALIGNMENTS COMPLETE"
fi
