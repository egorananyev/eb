# Trial exclusion

write_exclusions = function(trials_to_exclude, this_data_dir){
    # Logging which trials are marked for exclusion:
    write.table(trials_to_exclude, paste(this_data_dir, 'exclude_trials.csv', sep='/'),
                row.names=F, col.names=F, sep=',')
}

get_exclusions = function(this_data_dir){
    excluded_trials = read.csv(paste(this_data_dir, 'exclude_trials.csv', sep='/'),
                               header=F)
}
