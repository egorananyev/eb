# Reading ET functions
# Egor Ananyev
# 2019-01-12

# Getting the data directory
get_dir = function(dropbox_dir, eb, cue, subj, cond, sess){
    cond_dir = paste(paste0(dropbox_dir, 'Projects/eb/data/eb', eb), cue, subj, cond, sep='/')
    this_sess = paste0('block-', as.character(sess))
    data_dir = paste(cond_dir, dir(cond_dir, pattern=this_sess), sep='/')
    # TEMP: for now, just taking the last directory:
    # all_data_dirs = list.dirs(cond_dir)[-1]  # excluding the first (base) directory
    # data_dir = all_data_dirs[length(all_data_dirs)]
    return(data_dir)
}

# Function for reading faw .asc.gz files
read_em = function(data_dir){
    file_name = dir(data_dir, pattern='.asc.gz')
    print(file_name)
    print(dir(data_dir))
    file_path = paste(data_dir, file_name, sep='/')
    # Read in the data:
    print(file_path)
    raw_file = gzfile(file_path)
    raw_data = readLines(raw_file)  # gzfile unpacks .gz
    close(raw_file)
    return(raw_data)
}

# This function converts a specific chunk of raw data into a numeric data frame.
dfy = function(df_in, cols, space_col=0){
    ## The 'cols' parameter identifies the number of columns that are tab-delimited
    ## The 'space_col' parameter allows specification of _one_ column to delimit by space
    if(space_col == 0){
        # split the dataset by the tab '\t' delimeter:
        df_out = do.call(rbind, strsplit(df_in, '\t'))[,cols]
    }else{
        # split the dataset by the tab '\t' delimeter:
        df_out = data.frame(do.call(rbind, strsplit(df_in, "\t", fixed=TRUE)))
        # split the data set by space ' ' delimeter:
        df_out = cbind(df_out, data.frame(do.call(rbind, 
                                                  strsplit(as.character(df_out[,space_col])," ", 
                                                           fixed=TRUE))))
        df_out = df_out[,cols]
    }
    # We then convert the resulting character matrix into a numeric data frame:
    if(is.null(dim(df_out))){  # if there is only one dimension (i.e. list), no need f/ 'apply'
        df_out = as.numeric(df_out)
    }else{  # otherwise, 'apply'ing 'as.numeric' conversion to columns:
        if(space_col == 0){
            # To avoid the 'conversion to NA' warning, preassigning blink-related ROWS to NAs:
            df_out[df_out=='   .'] = NA
            df_out = data.frame(apply(df_out[,cols], 2, as.numeric))
        }else{
            df_out = data.frame(apply(df_out, 2, as.numeric))
        }
    }
    return(df_out)
}

parse_samples = function(raw_data){
    # extracting lines that start with numbers:
    samples_raw = raw_data[grepl('^[0-9]', raw_data)]
    # creating a proper data frame:
    samples = dfy(samples_raw, 1:4, 0)  # we only need 4 cols, and no columns are space-delimited
    colnames(samples) = c('sample', 'xr', 'yr', 'psr')
    print(head(samples))
    return(samples)
}

# Reading all once-per-trial events, including cue, blink latency, blink window, & response:
parse_trials = function(raw_data, cond){
    # Storing the start and end trials times:
    trials = data.frame(start=dfy(raw_data[grepl('TRIAL_START', raw_data)], c(3,5), 2),
                        end=dfy(raw_data[grepl('TRIAL_END', raw_data)], c(3,5), 2))
    trials$trial = as.numeric(rownames(trials))
    # Adding cue onset time:
    trials = cbind(trials, 
                   data.frame(cue_onset=dfy(raw_data[grepl('CUE_ONSET', raw_data)], c(3,5), 2)))
    # Adding blink latency onset time:
    trials = cbind(trials, 
                   data.frame(blink_latency=dfy(raw_data[grepl('BLINK_LATENCY_ONSET', raw_data)],
                                                c(3,5), 2)))
    # Adding blink window onset time:
    trials = cbind(trials, 
                   data.frame(blink_window=dfy(raw_data[grepl('BLINK_WINDOW_ONSET', raw_data)],
                                                c(3,5), 2)))
    # Adding target onset (confusingly called RESPONSE_ONSET in the eye-tracking output --
    # ... too late to change that, I suppose):
    trials = cbind(trials, 
                   data.frame(blink_window=dfy(raw_data[grepl('RESPONSE_ONSET', raw_data)], 
                                                c(3,5), 2)))
    if(cond!='cond-m'){
        # Adding response onset time:
        trials = cbind(trials, 
                       data.frame(blink_window=dfy(raw_data[grepl('TRIAL_RESPONSE', raw_data)], 
                                                    c(3,5), 2)))
        # Renaming columns to prettier variable names:
        colnames(trials) = c('trial_sample_beg', 'trial_time_beg', 
                             'trial_sample_end', 'trial_time_end', 'trial',
                             'cue_sample', 'cue_time', 'blink_latency_sample', 'blink_latency_time',
                             'blink_window_sample', 'blink_window_time', 'targ_sample', 'targ_time',
                             'resp_sample', 'resp_time')
    } else {
        # Renaming columns to prettier variable names:
        colnames(trials) = c('trial_sample_beg', 'trial_time_beg', 
                             'trial_sample_end', 'trial_time_end', 'trial',
                             'cue_sample', 'cue_time', 'blink_latency_sample', 'blink_latency_time',
                             'blink_window_sample', 'blink_window_time', 'targ_sample', 'targ_time')
        
    }
    # Adding the on/off times for the shutter goggles:
    if(cond == 'cond-a'){
        shutters_on = dfy(raw_data[grepl('SHUTTER_START', raw_data)], c(3,5), 2)
        # Due to a fixed bug, some shutter off fail to appear: accounting for that
        shutters_off = shutters_off=dfy(raw_data[grepl('SHUTTER_END', raw_data)], c(3,5), 2)
        # Registering the space bar presses in the AB condition:
        cue_resp_lines = grepl('CUE_RESP_TIME', raw_data)
        if(sum(cue_resp_lines) > 0){
            space_time = dfy(raw_data[cue_resp_lines], c(3,5), 2)
        } else {
            space_time = NA
        }
        # Labeling the shutter-off time trials:
        for(cur_trial in as.numeric(rownames(trials))){
            this_trial_sample_beg = trials$trial_sample_beg[cur_trial]
            this_trial_sample_end = trials$trial_sample_end[cur_trial]
            shutters_off$trial[shutters_off$X1 >= this_trial_sample_beg &
                               shutters_off$X1 <= this_trial_sample_end] = cur_trial
            shutters_on$trial[shutters_on$X1 >= this_trial_sample_beg &
                              shutters_on$X1 <= this_trial_sample_end] = cur_trial
            if(!is.null(nrow(space_time))){
                space_time$trial[space_time$X1 >= this_trial_sample_beg &
                                 space_time$X1 <= this_trial_sample_end] = cur_trial
            }
        }
        trials = merge(trials, shutters_off, by='trial', all=T)
        trials = merge(trials, shutters_on, by='trial', all=T)
        if(!is.null(nrow(space_time))){
            trials = merge(trials, space_time, by='trial', all=T)
        } else {
            trials$na1 = NA
            trials$na2 = NA
        }
        colnames(trials)[(ncol(trials)-5):ncol(trials)] = c('shutter_sample_beg','shutter_time_beg',
                                                            'shutter_sample_end','shutter_time_end',
                                                            'cue_resp_sample','cue_resp_time')
    }
    # For sample rate check, the last column should be all 1000 (corresponding to Hz):
    trials$tot_trial_samples = trials$trial_sample_end - trials$trial_sample_beg
    trials$tot_trial_time = trials$trial_time_end - trials$trial_time_beg
    trials$samples_per_s = trials$tot_trial_samples / trials$tot_trial_time
    print(paste('number of trials is', as.character(nrow(trials))))
    print(head(trials))
    return(trials)
}

# Labeling samples with trial numbers and adding 'time' column (in seconds)
lab_samples = function(samples, trials){
    all_samples = data.frame()
    for(cur_trial in as.numeric(rownames(trials))){
        this_trial_sample_beg = trials$trial_sample_beg[cur_trial]
        this_trial_sample_end = trials$trial_sample_end[cur_trial]
        trial_samples = samples[samples$sample>=this_trial_sample_beg & 
                                samples$sample< this_trial_sample_end,]
        trial_samples$trial = cur_trial
        trial_samples$time = (trial_samples$sample - this_trial_sample_beg) * 0.001
        all_samples = rbind(all_samples, trial_samples)
    }
    samples = all_samples
    print(head(samples))
    return(samples)
}

# Read blank events (may be more than one per trial)
parse_blanks = function(raw_data, trials){
    # Extracting the lines with blank ends, as they include blank start, end, and duration:
    blanks = dfy(raw_data[grepl('^EBLINK', raw_data)], c(6,2,3), 1)
    if(ncol(blanks)==1){
        blanks = t(blanks)
        rownames(blanks) = NULL
    }
    colnames(blanks) = c('blank_sample_beg', 'blank_sample_end', 'tot_blank_samples')
    all_blanks = data.frame()
    # Labeling trials in the 'blanks' data frame
    for(cur_trial in as.numeric(rownames(trials))){
        this_trial_sample_beg = trials$trial_sample_beg[cur_trial]
        this_trial_sample_end = trials$trial_sample_end[cur_trial]
        if(nrow(blanks)>1){
            trial_blanks = blanks[blanks$blank_sample_beg>=this_trial_sample_beg &
                                  blanks$blank_sample_beg< this_trial_sample_end,]  # ... note that
            # ... this implies that the blank must have been initiated on this trial
        } else {
            if(nrow(blanks)==1){  # if there's only one blank, see if it matches the trial
                if(blanks[1] >= this_trial_sample_beg & blanks[2] < this_trial_sample_end){
                    trial_blanks = blanks
                } else {
                    trial_blanks = NULL
                }
            } else {
                trial_blanks = NULL
            }
        }
        # The following can only be run if the trial_blanks is non-empty:
        if(!is.null(trial_blanks)){
            if(nrow(trial_blanks) > 0){
                trial_blanks = data.frame(trial_blanks, data.frame(trial = cur_trial))
                this_blank_time_beg = (trial_blanks$blank_sample_beg - this_trial_sample_beg) *  0.001
                trial_blanks = data.frame(trial_blanks, data.frame(blank_time_beg = this_blank_time_beg))
                this_blank_time_end = (trial_blanks$blank_sample_end - this_trial_sample_beg) * 0.001
                trial_blanks = data.frame(trial_blanks, data.frame(blank_time_end = this_blank_time_end))
                all_blanks = rbind(all_blanks, trial_blanks)
            }
        }
    }
    if(nrow(all_blanks) > 0){
        all_blanks$tot_blank_time = all_blanks$blank_time_end - all_blanks$blank_time_beg
    } else {
        all_blanks = data.frame(blank_sample_beg = NA, blank_sample_end = NA,
                                tot_blank_samples = NA, trial = NA,
                                blank_time_beg = NA, blank_time_end = NA)
    }
    print(head(all_blanks))
    return(all_blanks)
}
