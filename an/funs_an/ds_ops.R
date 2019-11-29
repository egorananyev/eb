## Egor Ananyev
## 2019-11-14

## ds_ops.R

## Data-set operations such as naming and variable scaling.

# Creating additional columns with level names that are easier to read:
ds_naming = function(ds){
  ds$Cond = 'NoBlink'
  ds$Cond[ds$cond=='a'] = 'Artificial'
  ds$Cond[ds$cond=='v'] = 'Prompted'
  ds$CondFull = 'No Blink'
  ds$CondFull[ds$cond=='a'] = 'Artificial Blink'
  ds$CondFull[ds$cond=='v'] = 'Prompted Blink'
  ds$CuePred = 'Predictive'
  ds$CuePred[ds$cue_pred==0] = 'Unpredictive'
  ds$CuePredFull = 'Predictive Cue'
  ds$CuePredFull[ds$cue_pred==0] = 'Unpredictive Cue'
  ds$CueValid = 'Valid'
  ds$CueValid[ds$cue_valid==0] = 'Invalid'
  return(ds)
}

# Rescaling and binning variables:
ds_vars = function(ds, ctoa_bin_nof){
  # Centered variables:
  ds$cue_valid_c = ds$cue_valid
  ds$cue_valid_c[ds$cue_valid==0] = -1
  ds$cue_pred_c = ds$cue_pred
  ds$cue_pred_c[ds$cue_pred==0] = -1
  # Reaction time:
  ds$RT = 1000 * ds$rt
  # CTOA:
  ds$ctoa_s = .1 + ds$blink_latency + .3 + ds$targ_soa
  ds$ctoa = ds$ctoa_s * 1000
  # Binning CTOA:
  ctoa_bins = seq(640, max(ds$ctoa), length.out=ctoa_bin_nof)  # min(ds$ctoa) -> 640
  ds$ctoa_bin = cut(ds$ctoa, ctoa_bins)
  ds$ctoa_lowbound = as.integer(ctoa_bins[as.numeric(ds$ctoa_bin)]) # lower bound
  print(paste('Excluding', as.character(sum(is.na(ds$ctoa_lowbound))), 'NA values due to binning.'))
  ds = ds[!is.na(ds$ctoa_lowbound),]  # excluding NAs
  ds$ctoa_lowbound_s = (ds$ctoa_lowbound / 1000) - (min(ds$ctoa_lowbound) / 1000)
  # We need a CTOA variable where CTOA_0[640] = 0
  ds$ctoa0_s = ds$ctoa_s - min(ds$ctoa_s)
  # Also adding the CTOA variable with three levels:
  ds$ctoa3 = 'Medium CTOA'
  ds$ctoa3[ds$ctoa0_s < max(ds$ctoa0_s)/3] = 'Short CTOA'
  ds$ctoa3[ds$ctoa0_s > 2 * max(ds$ctoa0_s)/3] = 'Long CTOA'
  # The 0-0.5-1 version of the above:
  ds$ctoa3_bin0 = 0
  ds$ctoa3_bin0[ds$ctoa3=='Medium CTOA'] = 0.5
  ds$ctoa3_bin0[ds$ctoa3=='Long CTOA'] = 1
  return(ds)
}

# Merging the behavioral data with the eye tracking data:
ds_et_merge = function(et, ds){
  # For analyzing eye blink parameters alongside behavioral data, only prompted blink cond is meaningful:
  et_filt = et[et$cond!='c',]
  # (This could have been done at the input stage, but I may need post-trial blinks at some point.)
  
  # In the <eb> dataset, only keeping blinks that occurs between the cue and the target:
  et_filt2 = et_filt[(et_filt$blank_sample_beg > et_filt$cue_sample) & 
                 (et_filt$blank_sample_end < et_filt$targ_sample) &
                 (!is.na(et_filt$blank_sample_beg)),]
  # CHECK: Test row to check if the above works:
  # ( et_filt[et_filt$subj==3 & et_filt$sess==1 & et_filt$cond=='v' & et_filt$cue_pred==1 & et_filt$trial==11,] )
  # ( et_filt2[et_filt2$subj==3 & et_filt2$sess==1 & et_filt2$cond=='v' & et_filt2$cue_pred==1 & et_filt2$trial==11,] )
  
  # Also excluding trials with more than one blink in the above time frame:
  num_blanks = ddply(et_filt2, .(subj, cue_pred, cond, sess, trial), 
                     summarise, num_blanks = length(trial))
  et_filt3 = anti_join(et_filt2, num_blanks[num_blanks$num_blanks>1,],
                      by=c('subj', 'cue_pred', 'cond', 'sess', 'trial'))
  # CHECK: Based on the above filter, the following number of blanks (rows) should have been removed:
  # sum(num_blanks[num_blanks$num_blanks>1,'num_blanks'])  # should be the same as:
  # nrow(et_filt2) - nrow(et_filt3)
  
  # Merging:
  ds = merge(ds[ds$cond!='c',], et_filt3, by=c('subj', 'cue_pred', 'sess', 'trial', 'cond'))
  return(ds)
}

ds_exclude = function(excl_ds, excl_reason_text){
  if(nrow(excl_ds) > 0){  # if there are any trials to exclude based on the short RTs
    excl_ds$excl_reason = excl_reason_text
    excl = rbind(excl, excl_ds)
    rownames(excl) <- NULL
    print(paste(excl_reason_text, as.character(nrow(excl_ds))))
  }
  return(excl)
}