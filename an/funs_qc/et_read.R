# Reading ET functions
# Egor Ananyev
# 2019-01-12

# Getting the data directory
get_dir = function(cloud_dir, eb, cue, subj, cond, sess){
    cond_dir = paste(paste0(cloud_dir, 'Projects/eb/data/eb', eb), cue, subj, cond, sep='/')
    if(!eb=='3'){
        this_sess = paste0('sess-', as.character(sess))
    } else {
        this_sess = paste0('block-', as.character(sess))
    }
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
    ## This identifies the blanks and puts them into a data frame along with start/end/duration,
    ## as well as whether the blank (blink) occurred post-cue and post-target.
    # Extracting the lines with blank ends, as they include blank start, end, and duration:
    eblink_indices = grepl('^EBLINK', raw_data)
    if(sum(eblink_indices) == 0){  # That would mean *no* blinks -- it's possible!
        all_blanks = data.frame(blank_sample_beg = NA, blank_sample_end = NA,
                                tot_blank_samples = NA, trial = NA,
                                blank_time_beg = NA, blank_time_end = NA, tot_blank_time = NA,
                                blank_post_cue = NA, blank_post_targ = NA)
        # all_blanks = NULL  # This doesn't work, as append fails for NULL vars
        print('No blinks detected: returning NA columns.')
    }else{
        blanks = dfy(raw_data[eblink_indices], c(6,2,3), 1)
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
                trial_blanks = blanks[blanks$blank_sample_beg >= this_trial_sample_beg &
                                      blanks$blank_sample_beg < this_trial_sample_end, ]
                # ... this implies that the blank must have been initiated on this trial
            } else {
                if(nrow(blanks)==1){  # if there's only one blank, see if it matches the trial
                    if(blanks[1] >= this_trial_sample_beg & blanks[2] < this_trial_sample_end){
                        trial_blanks = data.frame(blanks)
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
                    # trial_blanks = data.frame(trial_blanks, data.frame(trial = cur_trial))
                    trial_blanks$trial = cur_trial
                    # this_blank_time_beg = (trial_blanks$blank_sample_beg - this_trial_sample_beg) * 0.001
                    # trial_blanks = data.frame(trial_blanks,
                    #                           data.frame(blank_time_beg = this_blank_time_beg))
                    trial_blanks$blank_time_beg = 0.001 *
                        (trial_blanks$blank_sample_beg - this_trial_sample_beg)
                    # this_blank_time_end = (trial_blanks$blank_sample_end - this_trial_sample_beg) * 0.001
                    # trial_blanks = data.frame(trial_blanks,
                    #                           data.frame(blank_time_end = this_blank_time_end))
                    trial_blanks$blank_time_end = 0.001 *
                        (trial_blanks$blank_sample_end - this_trial_sample_beg)
                    # To see whether a given blank occurred post-cue and post-targ,
                    # ...taking the trial info from the <trials> df:
                    this_trial = trials[trials$trial==cur_trial, ]
                    trial_blanks$tot_blank_time = with(trial_blanks,
                                                       blank_time_end - blank_time_beg)
                    trial_blanks$blank_post_cue = as.numeric(trial_blanks$blank_sample_beg >
                                                       this_trial$cue_sample)
                    trial_blanks$blank_post_targ = as.numeric(trial_blanks$blank_sample_beg >
                                                        this_trial$targ_sample)
                    all_blanks = rbind(all_blanks, trial_blanks)
                }
            }
        }
        if(nrow(all_blanks) == 0){
            all_blanks = data.frame(blank_sample_beg = NA, blank_sample_end = NA,
                                    tot_blank_samples = NA, trial = NA,
                                    blank_time_beg = NA, blank_time_end = NA, tot_blank_time = NA,
                                    blank_post_cue = NA, blank_post_targ = NA)
            print('No blinks detected: returning NA columns.')
        } else {
            print(head(all_blanks))
        }
    }
    return(all_blanks)
}

# Read saccade events (may be more than one per trial) based on saccadic threshold:
label_sacc_samples = function(samples, blanks, trials, sacc_thresh = 2, eb_buff = 50, 
                              pre_cue_samples = 200){
    ### sacc_thresh = threshold for the HEM to exceed to count as saccade, in DOVA
    ### eb_buff = buffer time around an eye blink where saccades are ignored; in samples / ms
    
    # The only output from this is a modified <samples> variable, which contains the eye-tracking
    # data. The only new thing this function adds is the $sacc column that has a zero for
    # no saccade and one for saccadic sample.
    
    # First, setting the samples$sacc to all zeros:
    samples$sacc = 0
    # Also adding columns that will mark which trial phase the samples belong to:
    samples$pre_cue = NA  # recording hori eye post for pre_cue_samples prior to the cue
    samples$pre_eb = NA
    samples$post_eb = NA
    samples$pre_targ = NA  # this includes both pre_eb and post_eb in blink conditions
    samples$post_targ = NA
    samples$post_resp = NA
    
    # The following variable collects all samples labeled as saccades through the below loops:
    all_sacc_samples = NULL
    
    # cur_trial = 4  #TEMP
    for(cur_trial in trials$trial){
        # Taking the trial row from the <trials> df:
        this_trial = trials[trials$trial==cur_trial,]
        # All horizontal eye-tracking data for this trial:
        trial_etx = samples[samples$trial==cur_trial,]
        # Locating <pre_cue_samples> pre-cue samples:
        pre_cue_sam = trial_etx$sample[trial_etx$sample < this_trial$cue_sample &
                                       trial_etx$sample > (this_trial$cue_sample - pre_cue_samples)]
        pre_cue_ix = samples$sample %in% pre_cue_sam
        samples$pre_cue[pre_cue_ix] = samples$xr[pre_cue_ix]
        # Locating post-target (pre-response) samples:
        post_targ_sam = trial_etx$sample[trial_etx$sample > this_trial$targ_sample &
                                         trial_etx$sample < this_trial$resp_sample]
        post_targ_ix = samples$sample %in% post_targ_sam
        samples$post_targ[post_targ_ix] = samples$xr[post_targ_ix]
        # Locating post-response samples:
        post_resp_sam = trial_etx$sample[trial_etx$sample > this_trial$resp_sample]
        post_resp_ix = samples$sample %in% post_resp_sam
        samples$post_resp[post_resp_ix] = samples$xr[post_resp_ix]
        # For saccade search, limiting to eye-tracking samples after cue and before response:
        trial_etx = trial_etx[trial_etx$sample > trials$cue_sample[trials$trial==cur_trial],]
        trial_etx = trial_etx[trial_etx$sample < trials$resp_sample[trials$trial==cur_trial],]
        # We are only interested in pre-responce blanks for recording eye position:
        trial_blanks = blanks[blanks$trial == cur_trial &
                              blanks$blank_sample_end < this_trial$resp_sample, ]
        if(nrow(trial_blanks) == 0){
            no_blanks_found = TRUE
        } else {
            if(nrow(trial_blanks) == 1){
                if(is.na(trial_blanks[1,1])){
                    no_blanks_found = TRUE
                } else {
                    no_blanks_found = FALSE
                }
            } else {
                no_blanks_found = FALSE
            }
        }
        
        # Looping through blanks (if any) to label samples to examine -- only non-eye blink samples
        # (+/- a buffer zone, <buff_ms>) will be further examined:
        if(no_blanks_found == TRUE){  # if there are no blanks for this trial,
            trial_noneb_etx = trial_etx  # then all samples will be examined
        } else {  # if there are any blanks for this trial,
            # Since the following loops goes through however many blanks there are for this trial,
            # it is necessary to filter the same data set that many times:
            trial_noneb_etx = trial_etx
            pre_targ_eb_found = FALSE
            for(cur_blank in nrow(trial_blanks)){
                this_blank = trial_blanks[cur_blank,]
                pre_eb_inc = which(trial_noneb_etx$sample <
                                       (this_blank$blank_sample_beg - eb_buff))
                post_eb_inc = which(trial_noneb_etx$sample >
                                       (this_blank$blank_sample_end + eb_buff + 25))
                # If the blink happened before target appearance, recording the samples correspond-
                # ing to pre/post eb periods. It's only meaningful to do this once per trial, as
                # only one blink is allowed between cue and target (& only for some conds).
                if(!pre_targ_eb_found &
                   this_blank$blank_sample_beg > this_trial$cue_sample &
                   this_blank$blank_sample_end < this_trial$targ_sample) {
                    pre_eb_ix = samples$sample %in% trial_noneb_etx$sample[pre_eb_inc]
                    samples$pre_eb[pre_eb_ix] = samples$xr[pre_eb_ix]
                    post_eb_ix = samples$sample %in% trial_noneb_etx$sample[post_eb_inc]
                    samples$post_eb[post_eb_ix] = samples$xr[post_eb_ix]
                    # To label the pre-target samples, we still exclude eye blink samples:
                    samples$pre_targ[pre_eb_ix] = samples$xr[pre_eb_ix]
                    samples$pre_targ[post_eb_ix] = samples$xr[post_eb_ix]
                    # NA'ing the samples after the target appearance for $pre_targ:
                    samples$pre_targ[samples$sample >= this_trial$targ_sample] = NA
                    pre_targ_eb_found = TRUE
                }
                # Filtering the eye tracking data so that they only include pre/post-blank periods;
                # this is repeated for every blink:
                trial_noneb_etx = trial_noneb_etx[c(pre_eb_inc, post_eb_inc),]
            }
            # If there was no eye blink between cue and target, labeling the entire period as
            # pre-target:
            if(!pre_targ_eb_found){
                pre_targ_sam = trial_etx$sample[trial_etx$sample > this_trial$cue_sample &
                                                trial_etx$sample < this_trial$targ_sample]
                pre_targ_ix = samples$sample %in% pre_targ_sam
                samples$pre_targ[pre_targ_ix] = samples$xr[pre_targ_ix]
            }
        }
        
        # Recording the samples where the ET data exceed the saccadic threshold:
        trial_sacc_samples = trial_noneb_etx$sample[abs(trial_noneb_etx$xr) > sacc_thresh]
        
        # Appending the above to the saccade-collecting variable:
        all_sacc_samples = c(all_sacc_samples, trial_sacc_samples)
    }
    
    # Labeling the corresponding samples based on the above saccadic search:
    samples$sacc[samples$sample %in% all_sacc_samples] = 1
    
    return(samples)
}

# Measuring saccades for every trial.
# <sacc> = one saccade per row.
quant_sacc = function(samples, trials){
    # The following df is reduced in size every time a saccade is found for saccade search:
    subsamples = samples
    rownames(subsamples) <- NULL
    sacc = data.frame()
    
    # Making sure there are saccades in <samples>:
    perform_sacc_search = sum(samples$sacc==1) > 0  # TRUE or FALSE
    
    while(perform_sacc_search){
        # The row index of the beginning of the saccade:
        sacc_beg_ix = min(which(subsamples$sacc == 1))
        # The duration of the saccade, in number of rows (overshoots by 1):
        sacc_dur_rows = min(which(subsamples$sacc[sacc_beg_ix:nrow(subsamples)] == 0)) - 1
        # The row index of the end of the saccade:
        sacc_end_ix = sacc_beg_ix + sacc_dur_rows - 1
        
        this_sacc = data.frame(trial=NA,
                               sacc_beg_sample=NA, sacc_end_sample=NA, sacc_tot_sample=NA,
                               sacc_beg_time=NA, sacc_end_time=NA, sacc_tot_time=NA,
                               sacc_dir_R=NA, bef_targ=NA)
        this_trial = subsamples$trial[sacc_beg_ix]
        this_sacc$trial = this_trial
        # Recording the beginning, end, and duration of the saccade:
        this_sacc$sacc_beg_sample = subsamples$sample[sacc_beg_ix]
        this_sacc$sacc_end_sample = subsamples$sample[sacc_end_ix]
        this_sacc$sacc_tot_sample = with(this_sacc, sacc_end_sample - sacc_beg_sample)
        this_sacc$sacc_beg_time = subsamples$time[sacc_beg_ix]
        this_sacc$sacc_end_time = subsamples$time[sacc_end_ix]
        this_sacc$sacc_tot_time = with(this_sacc, sacc_end_time - sacc_beg_time)
        # A binary variable for whether the saccade is left (zero) or right (one):
        this_sacc$sacc_dir_R = as.numeric(subsamples$xr[sacc_beg_ix] > 0)
        # Classifying the saccade as occuring before or before target onset:
        this_sacc$bef_targ = as.numeric(this_sacc$sacc_beg_sample < 
                                            trials$targ_sample[trials$trial==this_trial])
        sacc = rbind(sacc, this_sacc)
        subsamples = subsamples[(sacc_end_ix+1):nrow(subsamples),]
        perform_sacc_search = sum(subsamples$sacc==1) > 0  # TRUE or FALSE
    }
    
    return(sacc)
}
