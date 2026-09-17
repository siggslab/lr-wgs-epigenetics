export OUTPUT_FILE=PGXXOX240048_pass_sup.sorted.retagged.nanomix.tsv
echo -e "chr\tstart\tend\ttotal_calls\tmodified_calls" > $OUTPUT_FILE
wgbstools beta2bed PGXXOX240048_pass_sup.sorted.retagged.beta | \
        awk -v OFS="\t" -v FS="\t" '{tmp=$5; $5=$4; $4=tmp; print}' >> $OUTPUT_FILE

