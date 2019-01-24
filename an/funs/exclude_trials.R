# Trial exclusion

write_exclusions = function(trials_to_exclude, data_dir){
    # Logging which trials are marked for exclusion:
    write.table(trials_to_exclude, paste(data_dir, 'exclude_trials.csv', sep='/'),
                row.names=F, sep=',')
}

# Manual override with additional trials to exclude (that weren't autodetected), or
# overriding falsely auto-excluded trials to be included back in.
get_manual = function(data_dir){
    # manual override -- exclusions (0) or inclusions (1);
    ## Example of a 'manual_override.csv' file content:
    # trial   excl0incl1    reason
    # 31      0 lost        signal
    ## Catalogue of possible options:
    # 0 "lost signal"
    # 0 "noisy"
    file_name = paste(data_dir, 'manual_override.csv', sep='/')
    if(file.exists(file_name)){
        manual = read.csv(paste(data_dir, 'manual_override.csv', sep='/'))
    } else {
        manual = NULL
    }
}
