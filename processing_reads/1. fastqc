# We use slurm and preambles for the scripts in these analyses are not provided. Those can be made available upon request

while read line || [ -n "$line" ]
do
        fastqc --noextract --nogroup --outdir path_to_fastq_files_dir $line
done < path_to_dir/my_file_paths.txt


#####################################################################################################################
#####################################################################################################################
############### Only need to change the path listed below. Nothing else needs to be changed to run this script.
###############
###############

path_to_files=""

###############
###############
###############
#####################################################################################################################
#####################################################################################################################


# Make new directory for the files if it doesn't exist (-p). Text files below will be overwritten if they exist.
mkdir -p fastQC_parse


echo -e "Position\tMean\tFile\tSample\tLane\tRead_number" > fastQC_parse/Per_base_sequence_quality.txt
echo -e "Position\tPercent\tBase\tFile\tSample\tLane\tRead_number" > fastQC_parse/Per_base_sequence_content.txt
echo -e "GC_Percent\tFragment_Count\tFile\tSample\tLane\tRead_number" > fastQC_parse/Per_sequence_GC_content.txt
echo -e "Sequence\tCount\tPercentage\tPossible_Source\tFile" > fastQC_parse/Overrepresented_sequences.txt
echo -e "Length\tCount\tFile\tSample\tLane\tRead_number" > fastQC_parse/Sequence_Length_Distribution.txt
echo -e "Duplication_Level\tPercent\tType\tFile\tSample\tLane\tRead_number" > fastQC_parse/Sequence_Duplication_Levels.txt


# Variable that controls the whole script. By default printing is turned off.
print="No"

# Loop through all fastQC zip files in the bear HiSeq folder. Variable i holds the path to each zip file each iteration through the loop.
# For each file, unzip the raw data (fastqc_data.txt) only and then write the relevant parts to separate files. The unzipped file
# gets written to the temp file (tteemmpp.crq) each iteration through the loop. Note that the file gets overwritten each time and then
# removed at the end. The if statement changes the value of the "print" variable. The case statement determines what information gets
# printed to what file depending on the value stored in the "print" variable.
for i in ${path_to_files}
do
        file_to_unzip=$(echo $i | awk -F "/" '{print $NF}')
        file_without_extension=${file_to_unzip%.*}
        # parse lane information, sample ID, and read number from filename
        # * in a regex is not like a filename glob. It means 0 or more of the previous character/pattern.
        # . in regex means "any character"
        lane_number=$(echo $file_to_unzip | grep -o "L00.")
        sample_ID=$(echo $file_to_unzip | awk -F "_" '{print $1}')
        read_number=$(echo $file_to_unzip | grep -o "R._" | sed 's/_//')


        unzip -p $i $file_without_extension"/fastqc_data.txt" > tteemmpp.crq


        while read line
        do
                if [[ ${line:2:10} == "END_MODULE" ]]; then
                        print="No"

                elif [[ ${line:2:25} == "Per base sequence quality" ]]; then
                        print="Per_base_sequence_quality"

                elif [[ ${line:2:25} == "Per base sequence content" ]]; then
                        print="Per_base_sequence_content"

                elif [[ ${line:2:23} == "Per sequence GC content" ]]; then
                        print="Per_sequence_GC_content"

                elif [[ ${line:2:25} == "Overrepresented sequences" ]]; then
                        print="Overrepresented_sequences"


                elif [[ ${line:2:28} == "Sequence Length Distribution" ]]; then
                        print="Sequence_Length_Distribution"

                elif [[ ${line:2:27} == "Sequence Duplication Levels" ]]; then
                        print="Sequence_Duplication_Levels"
                fi


                case $print in
                        "Per_base_sequence_quality")
                                echo "$line" | awk -v fn="$file_without_extension" \
                                                                   -v lane="$lane_number" \
                                                                   -v ID="$sample_ID" \
                                                                   -v read="$read_number" \
                                                                   '{print $1 "\t" $2 "\t" fn "\t" ID "\t" lane "\t" read}' >> fastQC_parse/Per_base_sequence_quality.txt
                                ;;

                        "Per_base_sequence_content")
                                echo "$line" | awk -v fn="$file_without_extension" \
                                                                   -v lane="$lane_number" \
                                                                   -v ID="$sample_ID" \
                                                                   -v read="$read_number" \
                                                                   '{print $1 "\t" $2 "\tG\t" fn "\t" ID "\t" lane "\t" read "\n"\
                                                                                   $1 "\t" $3 "\tA\t" fn "\t" ID "\t" lane "\t" read "\n"\
                                                                                   $1 "\t" $4 "\tT\t" fn "\t" ID "\t" lane "\t" read "\n"\
                                                                                   $1 "\t" $5 "\tC\t" fn "\t" ID "\t" lane "\t" read}' >> fastQC_parse/Per_base_sequence_content.txt
                                ;;

                        "Per_sequence_GC_content")
                                echo -e "$line\t$file_without_extension\t$sample_ID\t$lane_number\t$read_number" >> fastQC_parse/Per_sequence_GC_content.txt
                                ;;

                        "Overrepresented_sequences")
                                echo -e "$line\t$file_without_extension" >> fastQC_parse/Overrepresented_sequences.txt
                                ;;

                        "Sequence_Length_Distribution")
                                echo -e "$line\t$file_without_extension\t$sample_ID\t$lane_number\t$read_number" >> fastQC_parse/Sequence_Length_Distribution.txt
                                ;;

                        "Sequence_Duplication_Levels")
                                echo "$line" | awk -v fn="$file_without_extension" \
                                                                   -v lane="$lane_number" \
                                                                   -v ID="$sample_ID" \
                                                                   -v read="$read_number" \
                                                                   '{print $1 "\t" $2 "\tDeduplicated\t" fn "\t" ID "\t" lane "\t" read "\n"\
                                                                                   $1 "\t" $3 "\tTotal\t" fn "\t" ID "\t" lane "\t" read}' >> fastQC_parse/Sequence_Duplication_Levels.txt
                                ;;
                esac
        done < tteemmpp.crq
done

# Remove the temporary file used to parse.
rm tteemmpp.crq

# Remove the unwanted lines from the files (We added our own header above). GNU sed and masOS sed behave differently. You'll have to get rid
# of the '' behind -i for it to work on Kamiak. maOS: sed -i '' -e '/>>/d' -e '/#/d' >>>>>>> kamiak: sed -i -e '/>>/d' -e '/#/d'
sed -i -e '/>>/d; /#/d' fastQC_parse/Per_base_sequence_quality.txt
sed -i -e '/>>/d; /#/d' fastQC_parse/Per_base_sequence_content.txt
sed -i -e '/>>/d; /#/d' fastQC_parse/Per_sequence_GC_content.txt
sed -i -e '/>>/d; /#/d' fastQC_parse/Overrepresented_sequences.txt
sed -i -e '/>>/d; /#/d' fastQC_parse/Sequence_Length_Distribution.txt
sed -i -e '/>>/d; /#/d' fastQC_parse/Sequence_Duplication_Levels.txt


mv out_parse.txt fastQC_parse/
mv error_parse.txt fastQC_parse/
