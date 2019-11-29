# Blink estimation
# Egor Ananyev
# 2019-01-12

estimate_blinks = function(data_dir, exclusions){
    # To determine blink latency, merging blanks with beeps data sets:
    aggr = merge(trials, blanks, by='trial', all=T)
    # In measurement condition, excluding the blinks that occur after the trial ends:
    aggr = aggr[aggr$blank_time_end < (aggr$targ_time - aggr$trial_time_beg),]
    # Reading the trials to exclude
    if(exclusions){
        excluded_trials = read.csv(paste(data_dir, 'exclude_trials.csv', sep='/'))
        # Excluding the trials that didn't pass QC:
        aggr = aggr[!aggr$trial %in% excluded_trials, ]
    }
    # Computing blink latency:
    aggr$latency_samples = aggr$blank_sample_beg - aggr$cue_sample
    aggr$latency_time = aggr$blank_time_beg - aggr$cue_time + aggr$trial_time_beg
    # Computing blink duration:
    aggr$blink_duration = aggr$tot_blank_samples / 1000
    # Selecting columns for output and printing:
    select_cols = c('trial','latency_time','blink_duration')
    sub_df = aggr[,select_cols]
    write.csv(sub_df, paste(data_dir, 'blink_params.csv', sep='/'), row.names=F)
    print(sub_df)
}
