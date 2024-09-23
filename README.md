#  workflow for analysis of *Typha latifolia* in China

## Step01 : get the reference genome of *Typha latifolia*
here we get the reference from the study xxx
and the BioProject is PRJCA016992 in the link https://ngdc.cncb.ac.cn/gsub/submit/bioproject/list

After we get the reference genome, we next to do is the quality control of the raw data. The following is the file list in the raw data, and you also can check the data under the path data/china/china3/raw_data
    
    NGS data list: data.list
    

## Step02: quality control of the raw data
### Materials
    reference genome: Typha_latifolia.fasta (typha.fa)
    raw NGS data: see the data.list in data directory

### build the index of reference genome
    BWA: https://github.com/lh3/bwa
    Samtools: https://github.com/samtools/samtools
    GATK4: https://github.com/broadgsa/gatk
    bwa index -p typha typha.fa
    samtools faidx typha.fa
    java -Djava.io.tmpdir=/tmp -Xmx32g -jar path/gatk-package-4.2.4.1-local.jar CreateSequenceDictionary -R typha.fa -O typha.dict 

### quality control of the raw data
    fastp: https://github.com/OpenGene/fastp
    for sample in $(cat sample_list.txt);
    do 
    fastp -i ${sample}_1.fq.gz -I ${sample}_2.fq.gz -o ${sample}_1.fq.gz -O ${sample}_2.fq.gz -h ${sample}.html -j ${sample}.json -w 8 -q 20 -u 20 -n 5 -l 50 -y 20 -x --umi --umi_loc read1 --umi_len 10;
    done 

## Step03: mapping the clean data to the reference genome
    BWA: https://github.com/lh3/bwa
    Samtools: https://github.com/samtools/samtools
    for sample in $(cat sample_list.txt);
    do bwa mem -t 8 -R "@RG\tID:${sample}\tSM:${sample}\tPL:illumina" typha ${sample}_1.fq.gz ${sample}_2.fq.gz | samtools view -bS | samtools sort -o ${sample}.bam;
    samtools index ${sample}.bam;
    java -Djava.io.tmpdir=/tmp -Xmx128g -jar path/gatk-package-4.2.4.1-local.jar MarkDuplicates -I ${sample}.bam -O ${sample}_dedup.bam -M ${sample}_dedup.metrics.txt;
    samtools index ${sample}_dedup.bam;
    done

## Step04: call the SNP and InDel
    GATK4:  https://github.com/broadgsa/gatk
    for sample in $(cat sample_list.txt);
    do java -Djava.io.tmpdir=/tmp -Xmx128g -jar path/gatk-package-4.2.4.1-local.jar HaplotypeCaller --native-pair-hmm-threads 16 -R typha.fa -I ${sample}_dedup.bam -ERC GVCF -O ${sample}.g.vcf.gz;
    done

### joint genotyping into the all the samples
    GATK4:  https://github.com/broadgsa/gatk
    for file in $(cat list);do echo -en --variant ${file}.g.vcf' '; done > mergevcf
    echo  >> mergevcf
    cat mergevcf | while read line;do
    java -Djava.io.tmpdir=/tmp -Xmx128g -jar path/gatk-package-4.2.4.1-local.jar CombineGVCFs -R typha.fa $line -O merge.g.vcf.gz
    java -Djava.io.tmpdir=/tmp -Xmx128g -jar path/gatk-package-4.2.4.1-local.jar GenotypeGVCFs -R typha.fa -V merge.g.vcf.gz -O typha_raw.vcf.gz
    done
#### here we go the raw variants vcf file, and we need to do the filter of the variants
    typha_raw.vcf.gz
    

### joint genotyping of each collection site.


## Step05: filter the variants
    GATK4: https://github.com/broadgsa/gatk
    java -Djava.io.tmpdir=/tmp -Xmx128g -jar path/gatk-package-4.2.4.1-local.jar VariantsToTable -R typha.fa -V typha_raw.vcf.gz -O typha_raw_variants.csv -F CHROM -F POS -F DP -F AF -F QUAL -F QD -F MQ -F FS -F SOR -F MQRankSum -F ReadPosRankSum
    java -Djava.io.tmpdir=/tmp -Xmx128g -jar path/gatk-package-4.2.4.1-local.jar VariantFiltration -R typha.fa -V typha_raw.vcf.gz -O typha_raw_recode.vcf.gz --filter-expression "DP < 30||DP>8000 || QD < 15.0 || MQ < 58.0 || FS > 2.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -4.0" --filter-name "my_snp_filter"
    java -Djava.io.tmpdir=/tmp -Xmx128g -jar path/gatk-package-4.2.4.1-local.jar SelectVariants -R typha.fa -V typha_raw_recode.vcf.gz  -O typha_raw_filter_variants.vcf.gz --exclude-filtered true
    java -Djava.io.tmpdir=/tmp -Xmx128g -jar path/gatk-package-4.2.4.1-local.jar SelectVariants -R typha.fa -V typha_raw_filter_variants.vcf.gz --select-type-to-include SNP -O typha_raw_filter_variants_snp.vcf.gz

### here we got the final variants vcf file, and  we go the snp file 
    typha_raw_filter_variants.vcf.gz
    typha_raw_filter_variants_snp.vcf.gz

## Step06: build the phylogenetic tree of whole genome
    vcf2phylip: https://github.com/edgardomortiz/vcf2phylip
    IQ-TREE: https://github.com/Cibiv/IQ-TREE
    python /path/vcf2phylip.py -i typha_raw_filter_variants.vcf.gz -o typha_raw_filter_variants.phy
    iqtree -s typha_raw_filter_variants.phy -T 40 -bb 1000 -pre whole_genome_tree

#### after this step, we got the whole genome tree of the *Typha latifolia* in China and we visualize the tree in the iTOL
    iTOL: https://itol.embl.de/

## Step07: admixture analysis
    admixturePipeline: https://github.com/stevemussmann/admixturePipeline
    admixturePipeline.py -m popmap.txt -v typha_raw_filter_variants.vcf.gz -k 1 -K 10 -R 8  -n 48
#### here we use the CLUMPAK pipeline to basicly visualize the admixture results, and than use pytho script to plot the result
    python plot_admixture.py

## Step08: interface the history population size with PSMC
    PSMC: https://github.com/lh3/psmc
    for i in $(cat list)
    do
    samtools mpileup -C50 -uf typha.fa ${i}_markdup.bam | bcftools call -c | vcfutils.pl vcf2fq -d 10 -D 200 | gzip > ${i}.fq.gz
    fq2psmcfa -q20 ${i}.fq.gz > ${i}.psmcfa
    psmc -N25 -t15 -r5 -p '4+25*2+4+6' -o ${i}.psmc ${i}.psmcfa
    python script/plot_psmc.py -namelist popmap.txt -psmcdir psmc_data -o psmc_result.svg

## Step09: interface the history population size with momi2
    momi2: https://github.com/popgenmethods/momi2
#### here we prepare the python scrip for infering the history population size with momi2
    python script/momi_infer.py
#### and we check the log, to pick up the reasonable simulation result for the demographic history inference and we plot the resultw with python script
    python script/plot_momi.py

## Step10: build the whole chloroplast genome and the phylogenetic tree
    GetOrganelle: https://github.com/Kinggerm/GetOrganelle
    for i in `cat path_list`
    do 
    get_organelle_from_reads.py -1 $i_1.clean.fq.gz -2 $i_2.clean.fq.gz -o ./output/$i -R 15 -k 21,45,65,85 -F embplant_pt -t 8
    done
#### here we pick up a complete chloroplast genome as the reference to build the variant calling of the chloroplast genome
    BWA: https://github.com/lh3/bwa
    Samtools: https://github.com/samtools/samtools
    GATK4: https://github.com/broadgsa/gatk
    bwa index -p chloroplast chloroplast.fa
    samtools faidx chloroplast.fa
    java -Djava.io.tmpdir=/tmp -Xmx32g -jar path/gatk-package-4.2.4.1-local.jar CreateSequenceDictionary -R chloroplast.fa -O chloroplast.dict
    for sample in $(cat sample_list.txt);
    do bwa mem -t 8 -R "@RG\tID:${sample}\tSM:${sample}\tPL:illumina" chloroplast ${sample}_1.fq.gz ${sample}_2.fq.gz | samtools view -bS | samtools sort -o ${sample}.bam;
    samtools index ${sample}.bam;
    java -Djava.io.tmpdir=/tmp -Xmx128g -jar path/gatk-package-4.2.4.1-local.jar MarkDuplicates -I ${sample}.bam -O ${sample}_dedup.bam -M ${sample}_dedup.metrics.txt;
    samtools index ${sample}_dedup.bam;
    done
    
#### call the SNP and InDel

    for sample in $(cat sample_list.txt);
    do java -Djava.io.tmpdir=/tmp -Xmx128g -jar path/gatk-package-4.2.4.1-local.jar HaplotypeCaller --native-pair-hmm-threads 16 -R chloroplast.fa -I ${sample}_dedup.bam -ERC GVCF -O ${sample}.g.vcf.gz;
    done
### joint genotyping into the all the samples
    GATK4:  https://github.com/broadgsa/gatk
    for file in $(cat list);do echo -en --variant ${file}.g.vcf' '; done > mergevcf
    echo  >> mergevcf
    cat mergevcf | while read line;do
    java -Djava.io.tmpdir=/tmp -Xmx128g -jar path/gatk-package-4.2.4.1-local.jar CombineGVCFs -R chloroplast.fa $line -O merge.g.vcf.gz
    java -Djava.io.tmpdir=/tmp -Xmx128g -jar path/gatk-package-4.2.4.1-local.jar GenotypeGVCFs -R chloroplast.fa -V merge.g.vcf.gz -O typha_raw_chloroplast.vcf.gz
    done
#### here we got the raw variants vcf file, and we need to do the filter of the variants
    typha_raw_chloroplast.vcf.gz

### filter the variants
    GATK4:
    java -Djava.io.tmpdir=/tmp -Xmx128g -jar path/gatk-package-4.2.4.1-local.jar VariantsToTable -R chloroplast.fa -V typha_raw_chloroplast.vcf.gz -O typha_raw_chloroplast_variants.csv -F CHROM -F POS -F DP -F AF -F QUAL -F QD -F MQ -F FS -F SOR -F MQRankSum -F ReadPosRankSum
    java -Djava.io.tmpdir=/tmp -Xmx128g -jar path/gatk-package-4.2.4.1-local.jar VariantFiltration -R chloroplast.fa -V typha_raw_chloroplast.vcf.gz -O typha_raw_chloroplast_recode.vcf.gz --filter-expression "MQ < 40.0" --filter-name "MQ_filter" --filter-expression "MQRankSum < -1.5" --filter-name "MQRankSum_filter" --filter-expression "FS > 10.0" --filter-name "FS_filter" --filter-expression "QD < 15.0" --filter-name "QD_filter" --filter-expression "DP < 30.0" --filter-name "DP_large"
    java -Djava.io.tmpdir=/tmp -Xmx128g -jar path/gatk-package-4.2.4.1-local.jar SelectVariants -R chloroplast.fa -V typha_raw_chloroplast_recode.vcf.gz  -O typha_raw_chloroplast_filter_variants.vcf.gz --exclude-filtered true

#### here we got the final variants vcf file, and we build the phylogenetic tree of the chloroplast genome with IQ-TREE
    vcf2phylip:
    IQ-TREE:
    python /path/vcf2phylip.py -i typha_raw_chloroplast_filter_variants.vcf.gz -o typha_raw_chloroplast_filter_variants.phy
    iqtree -s typha_raw_chloroplast_filter_variants.phy -T 40 -bb 1000 -pre chloroplast_tree

#### after this step, we got the whole chloroplast genome tree of the *Typha latifolia* in China and we visualize the tree in the iTOL

###########################################################################

## Step11: sites' niche analysis
    Maxent: https://biodiversityinformatics.amnh.org/open_source/maxent/
    wordclim bio data: http://www.worldclim.com/past; https://www.worldclim.com/paleo-climate1;
    Last inter-glacial (LIG; ~120,000 - 140,000 years BP) : http://biogeo.ucdavis.edu/data/climate/worldclim/1_4/grid/pst/lig/lig_30s_bio.zip
    Last glacial maximum (LGM; ~21,000 years BP):http://biogeo.ucdavis.edu/data/climate/worldclim/1_4/grid/pst/21k/wc_2_5m_CCSM_21k_bio.zip
    Mid-Holocene (~6000 BP):http://biogeo.ucdavis.edu/data/climate/cmip5/mid/cnmidbi_2-5m.zip
    Current climate: https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_2.5m_bio.zip
    Future climate (2080-2100 ssp585): https://geodata.ucdavis.edu/cmip6/2.5m/MIROC6/ssp585/wc2.1_2.5m_bioc_MIROC6_ssp585_2081-2100.tif

### extract the environmental data of the sites and count the distribution and the principal component analysis of sites' environmental data
    python script/plot_pca_climate.py
    python script/plot_niche.py

### and the niche model of the *Typha latifolia* in China are build with 19 bioclimatic variables of above date and the Maxent software in GeoScence pro 

### and we visualize the niche model in the iTOL


## Step12: GWAS analysis
    vcf2gwas:https://github.com/frankvogt/vcf2gwas
    vcf2gwas -v typha_raw_filter_variants_snp.vcf.gz -pf env_data_current.csv -ap -lmm -T 6 -nl
    vcf2gwas -v typha_raw_filter_variants_snp.vcf.gz -pf env_data_6kya.csv -ap -lmm -T 6 -nl
    vcf2gwas -v typha_raw_filter_variants_snp.vcf.gz -pf env_data_21kya.csv -ap -lmm -T 6 -nl
#### as the position point of the data analysis of 130 kya are not fit into the on the map, so here we ignore the data analysis of 130 kya

### and we plot the GWAS result with the python script


## Step13: positive selection and selective sweep analysis
    RAiSD: https://github.com/alachins/raisd
    RAiSD -n typha -I typha_raw_filter_variants_snp.vcf.gz -O typha -R

### here we got the positive selection and selective sweep analysis result, and we plot the result with the python script


## Step14: RONA (risk of non-adaptive) analysis
#### here we use python script to calculate the RONA of the *Typha latifolia* in China


## Step15: gene offset analysis

