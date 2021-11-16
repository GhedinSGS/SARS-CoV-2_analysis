#! /bin/bash
#SBATCH -o /data/roderae/pipelines/covid/static/snpeffbuild%j.txt
#SBATCH -e /data/roderae/pipelines/covid/static/snpeffbuild-er%j.txt

set -e

module load snpEff
java -jar /usr/local/apps/snpEff/4.3t/snpEff.jar build -genbank -v SARS-CoV2_NC_045512.2 -c /data/roderae/pipelines/covid/static/snpEff.config
