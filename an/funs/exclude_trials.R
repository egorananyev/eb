# Trial exclusion

write_exclusions = function(trials_to_exclude){
    # Logging which trials are marked for exclusion:
    write.table(trials_to_exclude, paste(this_data_dir, 'exclude_trials.csv', sep='/'),
                row.names=F, col.names=F, sep=',')
}
