#To convert the VCF output to plink format
plink --vcf angsd2vcf_output_no_HG589.vcf.recode.vcf --recode --out plink_format/plink_formatangsd2vcf_output_no_HG589.vcf.recode --double-id --allow-extra-chr

#To get the bed bim fam file
plink --file plink_formatangsd2vcf_output_no_HG589.vcf.recode --allow-extra-chr --make-bed --out plink_formatangsd2vcf_output_no_HG589.vcf.recode_step1

#Edit the map file so that it has the correct format
sed -i 's/HiC_scaffold_//g' plink_formatangsd2vcf_output_no_HG589.vcf.recode.map

#Change column 2 to be the snp "name" in the map file
cat plink_formatangsd2vcf_output_no_HG589.vcf.recode.map | awk ' {print $1,"\t", $2=$1"_"$4,"\t",$3,"\t",$4} ' > plink_formatangsd2vcf_output_no_HG589.vcf.recode.map.edited

mv plink_formatangsd2vcf_output_no_HG589.vcf.recode.map.edited plink_formatangsd2vcf_output_no_HG589.vcf.recode.map

#Need to rename the ped file names

##First need to convert to a binary file format
plink --file plink_formatangsd2vcf_output_no_HG589.vcf.recode --make-bed --out mydata

cat mydata.fam | cut -f1 -d "_" #file is named rename.ids.txt

#Create a rename file (4 columns: oldFID oldIID newFID newIID)

Example (rename_ids.txt):

plink --bfile mydata --update-ids rename.ids.txt --make-bed --out mydata_renamed

#Convert back to the ped/map file format

plink --bfile mydata_renamed --recode --out mydata_renamed
