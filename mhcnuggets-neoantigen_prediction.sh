#! /usr/bin/env bash
#!/bin/bash

#abort on error
set -e

function usage
{
    echo "usage: neoantigen-prediction.sh -r <repo_root_path> -fq1 <fq1.fastq> -fq2 <fq2.fastq> -v <path.vcf> -t <threads> -o <outdir path>"
    echo "   ";
    echo "  -a | --hla_alleles          : If not genotyping: HLA genotypes in the format 'HLA-A02:01,HLA-B01:01'";
    echo "  -fq1 | --fastq_reads_1      : Fastq reads 1";
    echo "  -fq2 | --fastq_reads_2      : Fastq reads 2";
    echo "  -c | --class                : HLA class either I, II or both (Default both)";
    echo "  -v | --vcf                  : Path to VCF - NOTE cant be gzipped!";
    echo "  -r | --repo_root            : Path to repo root if not providing --arcasHLA and --mhcnuggets";
    echo "  --arcasHLA                  : Path to arcasHLA directory";
    echo "  --arcasHLA_env              : Name of conda env containing arcasHLA dependencies";
    echo "  --mhcnuggets                : Path to mhcnuggets directory";
    echo "  -t | --threads              : Number of threads (default 1)";
    echo "  -o | --outdir               : output directory path";

    echo "  -h | --help                 : This message";
}

function parse_args
{
  # positional args
  args=()

  # named args
  while [ "$1" != "" ]; do
      case "$1" in
          -r | --repo_root )            repo_root="$2";             shift;;
          -a | --hla_alleles )          hla_alleles="$2";     shift;;
          -fq1 | --fastq_reads_1 )      fq1="$2";     shift;;
          -fq2 | --fastq_reads_2 )      fq2="$2";      shift;;
          -c | --class )                class="$2";     shift;;
          -v | --vcf )                  vcf_path="$2";      shift;;
          -t | --threads )              threads="$2";      shift;;
          -o | --outdir )               outdir="$2";      shift;;
          --arcasHLA )                  arcasHLA="$2";      shift;;
          --mhcnuggets )                mhcnuggets="$2";      shift;;
          --arcasHLA_env )              arcasHLA_env="$2";      shift;;
          -h | --help )                 usage;                   exit;; # quit and show usage
          * )                           args+=("$1")             # if no match, add it to the positional args
      esac
      shift # move to next kv pair
  done

  # # restore positional args
  # set -- "${args[@]}"

  # # set positionals to vars
  # positional_1="${args[0]}"
  # positional_2="${args[1]}"

  # validate required args
  if [[ -z "$fq1" && -z "$fq1" && -z "$hla_alleles" ]]; then
    echo "Invalid arguments, need to provide --hla_alleles or -fq1 and -fq2"
    usage
    exit;
  fi

  if [[ -z "$outdir" || -z "$class" ]]; then
    echo "Missing required arguments"
    usage
    exit;
  fi

  if [[ -z "$arcasHLA" && -z "$repo_root" && -z "$mhcnuggets" ]]; then
    echo "Need to provide --repo_root or (--arcasHLA and --mhcnuggets)"
    usage
    exit;
  fi


  # set defaults
  if [[ "$class" == "I" ]]; then
    hla_genotype="A,B,C";
  elif [[ "$class" == "II" ]]; then
    hla_genotype="DMA,DMB,DOA,DOB,DPA1,DPB1,DQA1,DQB1,DRA,DRB1,DRB3,DRB5,E,F,G,H,J,K,L";
  elif [[ "$class" == "both" || -z "$class" ]]; then
    hla_genotype="A,B,C,DMA,DMB,DOA,DOB,DPA1,DPB1,DQA1,DQB1,DRA,DRB1,DRB3,DRB5,E,F,G,H,J,K,L";
  else
    echo "Invalid argument for --class"
    exit
  fi

  if [[ -z "$threads" ]]; then
    threads="1"
  fi

  if [[ -z "$arcasHLA" ]]; then
    arcasHLA="$repo_root/src/arcasHLA/arcasHLA"
  fi
  
  if [[ -z "$mhcnuggets" ]]; then
    mhcnuggets="$repo_root/src/mhcnuggets/mhcnuggets"
  fi

  if [[ -z "$arcasHLA_env" ]]; then
    arcasHLA_env="RA_analysis"
  fi

}



function summarize_errors
{
  #   Summarize errors in stdout
  grep "no mutation locat" $1/stdout.txt | cut -d"\(" -f1 > $1/stdout-summary-all.txt
  grep "alternative has more than one" $1/stdout.txt | cut -f1-7 -d" " >> $1/stdout-summary-all.txt
  grep "Invalid human allele" $1/stdout.txt -A1 | tr "\n" " " | sed 's/-- /\n/g' >> $1/stdout-summary-all.txt
  grep "ComplexSubstitution" $1/stdout.txt | cut -f1 -d"\(" >> $1/stdout-summary-all.txt
  
  echo -e "Count,Error-type" > $1/stdout-summary.csv
  sort < $1/stdout-summary-all.txt | uniq -c | sed -e 's/^ *//' | sed -r 's/\s+/,/' >> $1/stdout-summary.csv
  
  echo -e "Count,Error-type" > $1/stdout-summary-exclusions.csv
  grep -v -f $1/stdout-summary-all.txt -F $1/stdout.txt | sort | uniq -c | sed -e 's/^ *//' | sed -r 's/\s+/,/' >> $1/stdout-summary-exclusions.csv
  echo -e "Procesed errors in $1/stdout.txt - see $1/stdout-summary.csv for report"
  
}



function run
{
  parse_args "$@"

  mkdir -p $outdir
  
#  Run arcasHLA if needed
 if [[ ! -z "$fq1" && ! -z "$fq1" && -z "$hla_alleles" ]]; then
    >&2 echo -e "Activating conda environment: $arcasHLA_env..."
    eval "$(conda shell.bash hook)"
    conda activate $arcasHLA_env

    >&2 echo -e "\nRunning arcasHLA....\n"
    mkdir -p $outdir/arcas
    >&2 echo -e "\nRunning:\n$arcasHLA/arcasHLA genotype $fq1 $fq2 -g $hla_genotype -o $outdir/arcas -t $threads"
    $arcasHLA/arcasHLA genotype $fq1 $fq2 -g $hla_genotype -o $outdir/arcas -t $threads
    >&2 echo -e "\nRunning:\n$arcasHLA/arcasHLA merge --indir $outdir/arcas --outdir $outdir/arcas"
    $arcasHLA/arcasHLA merge --indir $outdir/arcas --outdir $outdir/arcas
    >&2 echo -e "\nRunnning:\n$arcasHLA/arcasHLA convert --resolution 2 $outdir/arcas/genotypes.tsv -o $outdir/arcas/genotypes-mod.tsv"
    $arcasHLA/arcasHLA convert --resolution 2 $outdir/arcas/genotypes.tsv -o $outdir/arcas/genotypes-mod.tsv
    >&2 echo -e "mv $outdir/arcas/genotypes-mod.tsv $outdir/hla_genotypes.tsv"
    mv $outdir/arcas/genotypes-mod.tsv $outdir/hla_genotypes.tsv

    hla_alleles=$(tail -n+2 < $outdir/hla_genotypes.tsv | tr "\\t" "," | sed "s/,/,HLA-/g" | cut -d"," -f2-15 | sed "s/\*//g")


  fi

  >&2 echo -e "\nActivating conda environment: mhcnuggets..."
  eval "$(conda shell.bash hook)"
  conda activate mhcnuggets
  



# Run class I prediction
if [[ "$class" == "I" || "$class" == "both" ]]; then
  
  # if runnning both cI and cII create subdir
  if [[ "$class" == "both" ]]; then
    outdir_run=$outdir/cI
  else 
    outdir_run=$outdir
  fi

  echo -e "\n\nRunning mhcnuggets.predict_from_vcf.py for class I - track stdout with tail -f $outdir_run/stdout.txt"

  mkdir -p $outdir_run

  # get cI alleles
  hla_geno=$(echo $hla_alleles | tr "," "\\n" | egrep "\-A|\-B|\-C" | sed -z 's/\n/,/g;s/,$/\n/')

  python $mhcnuggets/src/predict_from_vcf.py \
  -c I \
  -p I \
  -a $hla_geno \
  -g ensembl_grch38 \
  -v $vcf_path \
  -o $outdir_run > $outdir_run/stdout.txt 2> $outdir_run/stderr.txt

  summarize_errors $outdir_run

fi



# Run class II prediction
if [[ "$class" == "II" || "$class" == "both" ]]; then
  
  # if runnning both cI and cII create subdir
  if [[ "$class" == "both" ]]; then
    outdir_run=$outdir/cII
  else 
    outdir_run=$outdir
  fi

  echo -e "\n\nRunning mhcnuggets.predict_from_vcf.py for class II - track stdout with tail -f $outdir_run/stdout.txt"

  mkdir -p $outdir_run

  # get cII alleles
  hla_geno=$(echo $hla_alleles | tr "," "\\n" | egrep "\-DMA|\-DMB|\-DOA|\-DOB|\-DPA1|\-DPB1|\-DQA1|\-DQB1|\-DRA|\-DRB1|\-DRB3|\-DRB5|\-E|\-F|\-G|\-H|\-J|\-K|\-L" | sed -z 's/\n/,/g;s/,$/\n/')

  python $repo_root/src/mhcnuggets/mhcnuggets/src/predict_from_vcf.py \
  -c II \
  -p II \
  -a $hla_geno \
  -g ensembl_grch38 \
  -v $vcf_path \
  -o $outdir_run > $outdir_run/stdout.txt 2> $outdir_run/stderr.txt

  summarize_errors $outdir_run

fi


}

run "$@";


