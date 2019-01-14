# Blink estimation
# Egor Ananyev
# 2019-01-12

estimate_blinks = function(){
    # Tagging trial numbers in a separate column:
    df_beeps$trial = as.numeric(rownames(df_beeps))
    # To convert latency to seconds, we need to first compare sample count with measured time:
    df_beeps$tot_blink_samples = df_beeps$blank_sample_end - df_beeps$blank_sample_beg
    df_beeps$tot_time = df_beeps$blank_time_end - df_beeps$blank_time_beg
    df_beeps$sample_per_s = df_beeps$tot_time / df_beeps$tot_blink_samples
    # To detemine blink latency, merging blanks with beeps data sets:
    df_aggr = merge(blanks, df_beeps, by='trial')
    # Excluding the trials that didn't pass QC:
    df_aggr = df_aggr[!df_aggr$trial %in% trials_to_exclude, ]
    # Computing blink latency:
    df_aggr$latency_samples = df_aggr$blank_sample_beg.x - df_aggr$blank_sample_end.y
    df_aggr$latency_time = df_aggr$latency_samples * df_aggr$sample_per_s
    # Computing blink duration:
    df_aggr$blink_duration = df_aggr$tot_samples.x * df_aggr$sample_per_s
    # Selecting columns for output and printing:
    select_cols = c('trial','latency_time','blink_duration')
    sub_df = df_aggr[,select_cols]
    write.csv(sub_df, paste(this_data_dir, 'blink_params.csv', sep='/'), row.names=F)
    print(sub_df)
    return(sub_df)
}
