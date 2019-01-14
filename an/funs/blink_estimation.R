# Blink estimation
# Egor Ananyev
# 2019-01-12

estimate_blinks = function(this_data_dir){
    # Reading the trials to exclude
    excluded_trials = get_exclusions(this_data_dir)
    # To determine blink latency, merging blanks with beeps data sets:
    aggr = cbind(blanks, trials)
    # Excluding the trials that didn't pass QC:
    aggr = aggr[!aggr$trial %in% excluded_trials, ]
    # Computing blink latency:
    aggr$latency_samples = aggr$blank_sample_beg - aggr$cue_sample
    aggr$latency_time = aggr$blank_time_beg - aggr$cue_time + aggr$trial_time_beg
    # Computing blink duration:
    aggr$blink_duration = aggr$tot_blank_samples / 1000
    # Selecting columns for output and printing:
    select_cols = c('trial','latency_time','blink_duration')
    sub_df = aggr[,select_cols]
    write.csv(sub_df, paste(this_data_dir, 'blink_params.csv', sep='/'), row.names=F)
    print(sub_df)
}