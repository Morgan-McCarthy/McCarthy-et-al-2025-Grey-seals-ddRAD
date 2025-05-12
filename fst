module load angsd

#This is the reference genome you used for mapping (use the absolute path)
ref=/projects/mjolnir1/people/nrb613/README_gray_seal_2022/reference.genome/repeat_masked/Halichoerus_grypus_HiC_RM.fasta
#dir is where we will find the file containing a list of the bam files to be used. In my case I ran this where the bam files were but had the bam.list (containing the files to use) in a different directory (the calc_Fst directory). So the script will find the file that contains the names of bams to be used for a population in the directory indicated, but by default angsd will look for the actual files in the directory we are currently working in.
dir=cleaned_reads
#This is the output directory that is sort of the stem for all the analysis. It should be a directory where you run all the Fst analysis and contains bam.list files that have the bam files per location. I recommend keeping the calc_Fst ending because the script may refer back to it later on
out_dir=calc_Fst
#These are the chromosomes/scaffolds you are interested in. I have just isoloated the autosomes
#For me this file looks like this (just for the first 5) These should match the chromosome/scaffold headers in your reference genome
#SUPER_1
#SUPER_2
#SUPER_3
#SUPER_4
#SUPER_5
chrs=autosomes.txt
sites=final_filtered_global
#Now we start the loop after we have defined some global variables

for pop in $(cat unique.locations.txt)
do
#Here num is grabbing the number of individuals in a population.
#So num is a variable that reads out the file (cat) then it counts the number of entries and cuts the output based on fields 6 and 7 usng a space delimeter " ". It's just counting the number individuals in a population so we can use half that number as a filter in the next step
num=`cat ${dir}/${pop}_Bam_files.txt | wc | cut -f6,7 -d " "`

#This halfnum is a variable that takes the number of individuals in a file (num) and divides it by two. This way angsd will only consider sites that had data for at least 50% of the individuals within a population
halfnum=$(expr ${num} / 2)

#Now angsd will look for each ${pop}.bam.list and run the -dosaf 1 command with it as an input. 
#The -doSaf 1 command creates SAF files that are files that contain sample allele frequencies
#We also include the following filters/options 
#-anc ${ref} (this places the reference genome that we mapped everything to in the command. It's previously been saved as a variable before the loop)
#-gl 1 (use the samtools algorithm for genotype likelihoods)
#-minMapQ 30 (use a minimum mapping score of 30)
#-minQ 30 (use a minimum base quality score of 30)
#-minInd ${halfnum} (only consider sites where we had data for at least half of the samples in the bam list)
#-rf ${chrs} (only consider the chromosomes indicated here. I have only used the autosomes and excluded the sex chromosomes)
#-nThreads 8 (this is the number of threads)
#-out ${out_dir}${pop}/${pop} (this will put the .saf output in a directory named after the population and the prefix of the saf file will also be the name of the population. Note you don't need to add the suffix ".saf.idx".  
angsd -bam ${dir}/${pop}_Bam_files.txt -anc ${ref} -dosaf 1 -gl 1 -minMapQ 30 -minQ 30 -minInd ${halfnum} -rf ${chrs} -sites ${sites} -nThreads 8 -out ${out_dir}/${pop}/${pop}
done


##########################
Now that the input files are prepared we run
##########################

module load angsd 

for comparison in $(cat pop.compare.txt)
do
#Here I'm splling the two populations by their IDs using a ".". So pop1 will be a variable with the name of the first population and pop2 will be a variable with the name of the second population I'm interested in.

#These variables will be needed because for each population comparison ANGSD (and the loop) will write in files to directories it creates based on the populations

pop1="$(echo $comparison | cut -f1 -d"_")"
pop2="$(echo $comparison | cut -f2 -d"_")"

#This estimates the 2d-SFS
realSFS -P 10 ${pop1}/${pop1}.saf.idx ${pop2}/${pop2}.saf.idx > paired_2dsfs/$comparison.2dsfs

#This step creates an index file which will be needed
realSFS fst index -whichFst 1 ${pop1}/${pop1}.saf.idx ${pop2}/${pop2}.saf.idx -P 10 -sfs paired_2dsfs/${comparison}.2dsfs -fstout binary/${comparison}

#This step calculates the fst between the two populations in general (it gives you two values based on different calculations)
realSFS fst stats binary/${comparison}.fst.idx -P 20 > final_fsts/${comparison}.fst

#This step calculates the fst in nonoverlapping windows of 100 bps
realSFS fst stats2 binary/${comparison}.fst.idx -win 100 -step 100 -P 20 > final_fsts/${comparison}.100bp.windows.fst

done

