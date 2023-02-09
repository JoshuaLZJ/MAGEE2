PATH_TO_RELATE=/home/zhulaw/Documents/Honours/relate_v1.1.9_x86_64_static
vcf_file=/home/zhulaw/Documents/Honours/HGDP/hgdp_wgs.20190516.full.chrX.vcf.gz
#make sure the vcf.tbi file is also in the same directory where the vcf file is in!
position=75784694
chr=X
window_size=20000
pop=EAS #must be valid poplabel existing in samples_metadata.txt
genetic_map=/home/zhulaw/Documents/Honours/HGDP/genetic_map_chrX_nonPAR_combined_b37.txt

bcftools view -o windowed.vcf -r ${chr}:$(($position - $window_size / 2))-$(($position + $window_size / 2)) --threads 4 ${vcf_file} &> log_1.txt


mkdir ${pop}
python3 ./extract_samples.py -p ${pop}
python3 ./extract_poplabels.py -p ${pop}
cd ${pop}
bcftools view -S ${pop}_samples.txt ../windowed.vcf > ${pop}_windowed.vcf

bcftools filter -i 'TYPE="snp" && MIN(DP)>5 && QUAL>10' ${pop}_windowed.vcf > ${pop}_windowed_filtered.vcf
bcftools view --max-alleles 2 --exclude-types indels ${pop}_windowed_filtered.vcf > ${pop}_windowed_biallelic_filtered.vcf

bcftools convert ${pop}_windowed_biallelic_filtered.vcf --gensample ${pop}_windowed_biallelic_filtered

#genotype phasing
shapeit --input-gen ${pop}_windowed_biallelic_filtered.gen ${pop}_windowed_biallelic_filtered.sample -M ${genetic_map} -O ${pop}_windowed_biallelic_filtered.phased --chrX --force

bcftools convert ${pop}_windowed_biallelic_filtered.phased  --hapsample ${pop}_windowed

mv ${pop}_windowed.hap.gz ${pop}_windowed.haps.gz
gzip ${pop}_windowed.sample

${PATH_TO_RELATE}/bin/Relate --mode All --haps $PWD/${pop}_windowed.haps.gz --sample $PWD/${pop}_windowed.sample.gz --map ${genetic_map} -m 1.1e-8 -N 30000 -o ${pop}_windowed


${PATH_TO_RELATE}/scripts/TreeView/TreeViewMutation.sh --haps $PWD/${pop}_windowed.haps.gz --sample $PWD/${pop}_windowed.sample.gz --anc $PWD/${pop}_windowed.anc --mut $PWD/${pop}_windowed.mut --poplabels $PWD/${pop}.poplabels --bp_of_interest ${position} --years_per_gen 28 -o ${pop}_windowed

mkdir popsizes
cd popsizes

${PATH_TO_RELATE}/bin/RelateCoalescentRate --mode EstimatePopulationSize -i $PWD/${pop}_windowed -o ${pop}_windowed --bins 3,7,0.25 --poplabels ../${pop}.poplabels

cd ..

${PATH_TO_RELATE}/scripts/EstimatePopulationSize/EstimatePopulationSize.sh -i ${pop}_windowed -o ${pop}_windowed_ne --poplabels $PWD/${pop}.poplabels -m 1.1e-8 --years_per_gen 28 --bins 3,7,0.25 --num_iter 2 --threads 4

mkdir selection
cd selection

cp ../${pop}_windowed_ne* .

${PATH_TO_RELATE}/bin/RelateSelection --mode Frequency -i ${pop}_windowed_ne -o ${pop}_windowed_ne
${PATH_TO_RELATE}/bin/RelateSelection --mode Selection -i ${pop}_windowed_ne -o ${pop}_windowed_ne

cd ../..

export pop=$pop #for the relate.R script
Rscript ./relate.R