###!/usr/bin/arch -i386 /usr/bin/python # -*- coding: utf-8 -*-
# from __future__ import print_function

"""
The influence of eye blinks on attentional cueing.
Creator: Egor Ananyev
Original date: 2018-11-15

# Variables for communicating with Arduino (shutter goggles):
AirOff = 'a'
AirOn = 'b'
LeftOn = 'l' #goggles 
LeftOff = 'm' #goggles
RightOn = 'r' #goggles
RightOff = 's' #goggles
LensOn = 'c' #both sides clear
LensOff = 'z' #both sides opaque
AllOff = 'x'
"""

# from __future__ import division  # so that 1/3=0.333 instead of 1/3=0
from psychopy import prefs
prefs.hardware['audioLib'] = ['pyo']
# ['PTB'] fails after one run, but can be reset...
# ['sounddevice'] doesn't work
# ['pyo'] works, but with lots of ALSA warnings
import psychtoolbox as ptb
from psychopy import visual, core, event, gui, monitors, sound
import numpy as np
import os
from subprocess import call  # for running shell commands from within Python
from psychopy.data import TrialHandler, importConditions
from psychopy.core import wait
import pandas as pd
from datetime import datetime
import math

# Imports associated with the eye tracker:
import pylink
from EyeLinkCoreGraphicsPsychoPy import EyeLinkCoreGraphicsPsychoPy

## Initial variables.
# experiment modes:
toshi = False
dummy_mode = False
debug = True
drift_check = False
# experiment variables:
cue_delay_min = 500  # the time frame before the location/blink cue
cue_delay_max = 1000  # shortened from 800 to 500 ms on 2019-06-11
# display dimensions:
if toshi:
    # dr = (576, 432)
    ds = 61
    # dd = (17.2, 12.9)  # display dimensions in cm
    dr = (1152, 864)  # display resolution in px
    dd = (34.4, 25.8)  # display dimensions in cm
else:
    dr = (1152, 864)  # display resolution in px
    ds = 65  # distance to screen in cm
    dd = (40.0, 30.0)  # display dimensions in cm: 39.0x29.5; ~34deg horizontally, 1cm~0.85deg, 1deg~33.87px
# fixation cross:
fix_size = 0.8
background_color = [-.5, -.5, -.5]
# cue:
cue_size = 5
cue_dur = .1  # shortened on 2019-06-11
scue_dur = cue_dur  # will keep the spatial cue duration same as cue_dur for <eb4>
cue_color = [-1, -1, -1]
# target:
targ_off_x = 8
targ_diam = .8
targ_color = [0, 0, 0]

## getting user info about the experiment session:
exp_info = {u'subj': u'1', u'cond': u'c', u'block': u'', u'soa': u'00'}
# conditions: 't'=training, 'c'=control, 'a'=artificial blink, 'v'=voluntary blink, 'm'=measurement
# cue_pred: cue is either predictive (75% valid) or unpredictive (50% valid)
dlg = gui.DlgFromDict(dictionary=exp_info, title='eb')  # dialogue box
exp_info['expt'] = '4'
if not dlg.OK:
    core.quit()  # user pressed cancel
exp_info['time'] = datetime.now().strftime('%Y-%m-%d_%H%M')
exp_name = 'eb' + exp_info['expt']
print('experiment name is ' + exp_name)

cue_pred = None
trial_n = None
soa = int(exp_info['soa'])
if exp_name == 'eb1':
    # Predictiveness of the cue:
    cue_pred = int(exp_info['cue_pred'])
    trial_n = 5  # trials per condition row; 5 gives 40 trials (per block)
elif exp_name == 'eb2':
    cue_pred = 1
    trial_n = 1  # due to many SOA levels, only a single iteration can be performed per condition per block
elif exp_name == 'eb3':
    cue_pred = 1
    trial_n = 5  # trials per condition row; 5 gives 40 trials (per block); 20 ~ 9 min
    # with 11 blocks, the above session should last 40-50 min.
elif exp_name == 'eb4':
    cue_pred = 1
    trial_n = 2  # 2 trials per 24 conditions = 48 trials per block; 3 => 72

# Assigning conditions:
eye_tracking = True  # true by default
training = False
voluntary = False
shutters = False
measure = False
noblink = False
noblink_spacebar_recorded = False
print('Condition: ' + exp_info['cond'])
if exp_info['cond'] == 'm':
    measure = True
else:
    if exp_info['cond'] == 't':
        trial_n = 1
        training = True
    if exp_info['cond'] == 'v':
        voluntary = True
    if exp_info['cond'] == 'a':
        import serial
        shutters = True
        if not dummy_mode:
            ser = serial.Serial('/dev/ttyACM0', 9600)
            # ser.write('c')  # both sides clear
            ser.write(str.encode('c'))
        shutters_shut = False  # this keeps track of whether the shutters were shut
        shutters_opened = False  # this is an additional tracker for when the shutters were already closed on this trial
    if exp_info['cond'] == 'c':
        noblink = True
        shutters_shut = False
        shutters_opened = False

# Handling condition instructions:
if measure:
    cond_instr = 'Please:\n' \
                 '- Keep your head still;\n' \
                 '- Always look directly at the center cross (+);\n' \
                 '- When you see the BLACK SQUARE, blink once.'
elif exp_info['cond'] == 'c' or exp_info['cond'] == 'a':
    cond_instr = 'Please do the following:\n' \
                 '- Keep your head still and always look directly at the center cross;\n' \
                 '- DO NOT BLINK during the trial;\n' \
                 '- When you see the BLACK SQUARE, press ''SPACE'' immediately;\n' \
                 '- Pay attention to the beep:\n' \
                 '  * low-pitch = likely left target\n' \
                 '  * high-pitch = likely right target;\n' \
                 '- When you see the WHITE DOT, press ''LEFT'' or ''RIGHT'' to indicate its location.'
else:  # for 'd' or 'v'
    cond_instr = 'Please do the following:\n' \
                 '- Keep your head still and always look directly at the center cross;\n' \
                 '- When you see the BLACK SQUARE, BLINK ONCE immediately;\n' \
                 '- Pay attention to the beep:\n' \
                 '  * low-pitch = likely left target\n' \
                 '  * high-pitch = likely right target;\n' \
                 '- When you see the WHITE DOT, press ''LEFT'' or ''RIGHT'' to indicate its location.'

## Input and output

# Condition file:
exp_conditions = None
if not measure:
    if exp_name == 'eb1':
        if exp_info['cond'] == 'd':
            exp_conditions = importConditions('cond-files/cond_' + exp_name + '_d' + '.xlsx')
        else:
            # same design for all non-d conditions:
            if cue_pred:
                exp_conditions = importConditions('cond-files/cond_' + exp_name + '_cue_predictive.xlsx')
            else:
                exp_conditions = importConditions('cond-files/cond_' + exp_name + '_cue_unpredictive.xlsx')
    elif exp_name == 'eb2':  # for experiments eb2 and 3, a single conditions file is taken:
        exp_conditions = importConditions('cond-files/cond_' + exp_name + '.xlsx')
    elif exp_name == 'eb3':  # for eb3, SOA is determined by user input
        exp_conditions = importConditions('cond-files/cond_eb1_cue_predictive.xlsx')
    elif exp_name == 'eb4':  # for eb4, CTOA is determined by the user, but the cue-to-cue delay is from the condfile
        exp_conditions = importConditions('cond-files/cond_eb4.xlsx')
else:
    exp_conditions = importConditions('cond-files/cond_eb1_m.xlsx')

# The output directory will depend on whether the cue is predictive (regardless of whether it's eb1 or eb2):
if cue_pred:
    cue_dir = 'cp'
else:
    cue_dir = 'cu'

# Trial handler depending on the measure or experimental stage:
if measure:
    out_file_name = 'measure'
    trials = TrialHandler(exp_conditions, 12, extraInfo=exp_info)
else:
    out_file_name = 'beh_out'
    trials = TrialHandler(exp_conditions, trial_n, extraInfo=exp_info)

# output file:
exp_dir = '..' + os.sep + 'data' + os.sep + exp_name + os.sep + cue_dir
if not os.path.exists(exp_dir):
    print('experiment directory does not exist')
    os.makedirs(exp_dir)
else:
    print('experiment directory exists')
subj_dir = exp_dir + os.sep + 'subj-%02d' % int(exp_info['subj'])
if not os.path.exists(subj_dir):
    os.makedirs(subj_dir)
cond_dir = subj_dir + os.sep + 'cond-' + exp_info['cond']
if not os.path.exists(cond_dir):
    os.makedirs(cond_dir)
sess_dir = cond_dir + os.sep + 'block-%s_%s' % (exp_info['block'], exp_info['time'])
if not os.path.exists(sess_dir):
    os.makedirs(sess_dir)
out_file_path = sess_dir + os.sep + out_file_name + '.csv'

# output matrix:
output_mat = {}

## Shutters condition also estimates blink duration to simulate physiological blinks
if shutters:
    # taking only the last measurement:
    all_sess_dirs = os.listdir(subj_dir + os.sep + 'cond-m')
    last_sess_dir = all_sess_dirs[len(all_sess_dirs) - 1]
    blink_params = pd.read_csv(subj_dir + os.sep + 'cond-m' + os.sep + last_sess_dir + os.sep + 'blink_params.csv')
    blink_dur_ave = np.mean(blink_params['blink_duration'])
    blink_dur_std = np.std(blink_params['blink_duration'])
    # The later functions will run as follows:
    # np.random.normal(blink_dur_ave, blink_dur_std)

## EyeLink setup
if not dummy_mode:
    tracker = pylink.EyeLink('100.1.1.1')
else:
    tracker = pylink.EyeLink(None)

# Note that the file name cannot exceeds 8 characters. Open eyelink data files as early to record as possible.
edf_data_file_name = 'eye_out.edf'
tracker.openDataFile(edf_data_file_name)
# add personalized data file header (preamble text)
tracker.sendCommand("add_file_preamble_text 'Study: Influence of blinks on attentional cueing'")

## Monitor setup
if toshi:
    mon = monitors.Monitor('Toshi', width=dd[0], distance=ds)
    mon.setSizePix(dr)
    window = visual.Window(dr, monitor=mon, fullscr=False, screen=0, color=background_color, units='deg')
else:
    # you MUST specify the physical properties of your monitor first, otherwise you won't be able to properly use
    # different screen "units" in psychopy. One may define his/her monitor object within the GUI, but
    # I find it is a better practice to put things all under control in the experimental script instead.
    mon = monitors.Monitor('station3')  # , width=dd[0], distance=ds)
    # mon.setSizePix(dr)
    # print('Monitor distance')
    # print(mon.getDistance())
    window = visual.Window(dr, fullscr=True, monitor=mon, color=background_color, units='deg',
                           allowStencil=True, autoLog=False, screen=0, waitBlanking=False)

if toshi:
    frame_rate = 60
else:
    frame_rate = window.getActualFrameRate()
    if frame_rate < 99:
        print('WARNING! The measured frame rate is lower than expected')
print('frame rate: ' + str(frame_rate))

## Initialize the stimuli and instructions
space_text = "\n\nPress the spacebar to start"
instr_text = cond_instr + space_text
instr_text_stim = visual.TextStim(window, text=instr_text, height=.8)
fix_cross = visual.TextStim(window, text='+', bold='True', pos=[0, 0], rgb=1, height=fix_size)

# A box (blink) cue:
bcue_box = visual.Rect(window, width=cue_size, height=cue_size, lineColor=cue_color)

# A pitch (spatial) cue:
sample_rate = 16000  # 22500
freq_L = 200
freq_R = 800
# bits = -16

n_samples = int(round(cue_dur * sample_rate))
buf_L = np.zeros((n_samples, 2))
buf_R = np.zeros((n_samples, 2))
max_sample = 1.0

for s in range(n_samples):
    t = float(s) / sample_rate  # time in seconds
    val_L = max_sample * math.sin(2 * math.pi * freq_L * t)  # spans from -1 to 1
    val_R = max_sample * math.sin(2 * math.pi * freq_R * t)  # spans from -1 to 1
    # val = max_sample * val
    buf_L[s][0] = val_L  # left
    buf_L[s][1] = 0  # right
    buf_R[s][0] = 0  # left
    buf_R[s][1] = val_R  # right

print(sound.Sound)
scue_lo = sound.Sound(value=buf_L, sampleRate=sample_rate)  # , bits=bits)
scue_hi = sound.Sound(value=buf_R, sampleRate=sample_rate)  # , bits=bits)
print('loaded sound')

# A circle target:
targ = visual.Circle(window, radius=targ_diam / 2, edges=32, pos=(0, 5), fillColor=targ_color, lineColor=targ_color)

if not dummy_mode:
    ## Eye-tracking calibration:
    # call the custom calibration routine "EyeLinkCoreGraphicsPsychopy.py", instead of the default
    # routines that were implemented in SDL
    custom_calibration = EyeLinkCoreGraphicsPsychoPy(tracker, window)
    pylink.openGraphicsEx(custom_calibration)

    ## STEP V: Set up the tracker
    # we need to put the tracker in offline mode before we change its configurations
    tracker.setOfflineMode()
    # sampling rate, 250, 500, 1000, or 2000; this command won't work for EyeLInk II/I
    tracker.sendCommand('sample_rate 500')
    # inform the tracker the resolution of the subject display
    # [see Eyelink Installation Guide, Section 8.4: Customizing Your PHYSICAL.INI Settings ]
    tracker.sendCommand("screen_pixel_coords = 0 0 %d %d" % (dr[0] - 1, dr[1] - 1))
    # save display resolution in EDF data file for Data Viewer integration purposes
    # [see Data Viewer User Manual, Section 7: Protocol for EyeLink Data to Viewer Integration]
    tracker.sendMessage("DISPLAY_COORDS = 0 0 %d %d" % (dr[0] - 1, dr[1] - 1))
    # specify the calibration type, H3, HV3, HV5, HV13 (HV = horizontal/vertical),
    tracker.sendCommand("calibration_type = HV5")  # tracker.setCalibrationType('HV9') also works, see the Pylink manual
    # Set the tracker to parse Events using "GAZE" (or "HREF") data
    tracker.sendCommand("recording_parse_type = GAZE")
    # Online parser configuration: 0-> standard/cognitive, 1-> sensitive/psychophysiological
    # the Parser for EyeLink I is more conservative, see below
    # [see Eyelink User Manual, Section 4.3: EyeLink Parser Configuration]
    tracker.sendCommand('select_parser_configuration 0')
    # Host tracking software version is 5

    ## Sending eye-tracker commands and calibration:
    # specify the EVENT and SAMPLE data that are stored in EDF or retrievable from the Link
    # See Section 4 Data Files of the EyeLink user manual
    tracker.sendCommand("file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT")
    tracker.sendCommand("link_event_filter = LEFT,RIGHT,FIXATION,FIXUPDATE,SACCADE,BLINK,BUTTON,INPUT")
    tracker.sendCommand("file_sample_data  = LEFT,RIGHT,GAZE,AREA,GAZERES,STATUS,HTARGET,INPUT")
    tracker.sendCommand("link_sample_data  = LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS,HTARGET,INPUT")

    # Calibration:
    tracker.doTrackerSetup()

## Handy routines:

# Frame-skipping check:
def frame_skip_check(elapsed_t, elapsed_frames):
    # The number of elapsed frames should match the time:
    print('time=%.3f  frames=%d  rate=%.4f' % (elapsed_t, elapsed_frames, (elapsed_t / elapsed_frames)))

# This is done at every frame update, regardless of trial phase, so predefining:
def frame_routine():
    flip_time_ = window.flip()
    fix_cross.draw()
    # Checking for quit (the Esc key)
    if event.getKeys(keyList=['escape']):
        exit_routine()
    return flip_time_

## Artificial blink condition:
# Monitoring the 'spacebar' press in artificial and no-blink conditions:
def monitor_cue_resp(flip_time_, cue_rt_start_):
    space_key = event.getKeys(keyList=['space'])
    if len(space_key) > 0:
        print('The "space" key was pressed.')
        cue_rt_ = flip_time_ - cue_rt_start_
        return cue_rt_
    else:
        return 0

# Shutter-closing routine
def close_shutters(dummy_mode_):
    if not dummy_mode_:
        tracker.sendMessage('SHUTTER_START %.2f' % flip_time)
        # ser.write('z')
        ser.write(str.encode('z'))
    print('Closed the goggles.')
    shutters_shut_ = True
    return shutters_shut_

# Shutter-opening routine:
def open_shutters(dummy_mode_):
    if not dummy_mode_:
        tracker.sendMessage('SHUTTER_END %.2f' % flip_time)
        # ser.write('c')
        ser.write(str.encode('c'))
    print('Opened the goggles.')
    shutters_shut_ = False
    return shutters_shut_

# Also no variation across frames, but only available upon call, which is made only in key registering phase.
def exit_routine():
    # attempting to close the goggles:
    print('Initiating the exit routine.')
    # noinspection PyBroadException
    if shutters:
        try:
            # ser.write('z')
            ser.write(str.encode('z'))
            print('Closed the goggles.')
        except:
            print('Failed to close the goggles -- they might not have been used.')

    # Behavioural data output:
    if not output_mat:  # means that the dictionary is empty
        print('The output file is empty.')
    else:
        if measure:
            data_columns = ['exp_name', 'subj', 'cond', 'block', 'trial_id', 'cue_delay', 'trial_start', 'trial_end']
        else:
            data_columns = ['exp_name', 'cue_pred', 'subj', 'cond', 'block', 'trial_id', 'cue_delay', 'targ_right',
                            'scue_delay', 'cue_valid', 'targ_soa',
                            'shutter_dur', 'trial_start', 'trial_end', 'cue_rt', 'corr_resp', 'rt']
        pd.DataFrame.from_dict(output_mat, orient='index').to_csv(out_file_path, index=False, columns=data_columns)
        print('output file path is ' + out_file_path)

    # Say goodbye:
    window.flip()
    instr_text_stim.setText('    Finished!\nRecording data...')
    instr_text_stim.draw()
    window.flip()

    if not dummy_mode:
        # In case the EyeLink recording is ongoing, stopping the recording:
        pylink.pumpDelay(100)
        tracker.stopRecording()
        # EDF output:
        # close the EDF data file
        tracker.setOfflineMode()
        tracker.closeDataFile()
        pylink.pumpDelay(50)
        # Get the EDF data
        tracker.receiveDataFile(edf_data_file_name, sess_dir + os.sep + edf_data_file_name)
        # Converting the EDF to ASC from within this code
        print('converting EDF to ASC, zipping it, and moving the original EDF file to a different directory')
        call(['edf2asc', sess_dir + os.sep + edf_data_file_name])
        call(['gzip', sess_dir + os.sep + 'eye_out.asc'])
        call(['mv', sess_dir + os.sep + edf_data_file_name,
              '..' + os.sep + 'edf_data' + os.sep + exp_info['time'] + '.edf'])
        # Close the link to the tracker
        tracker.close()
        # Close the graphics
        pylink.closeGraphics()
    window.close()
    core.quit()


## Initiating the trial loop
n_trials_done = 0

for trial in trials:

    ## First trial initiates instructions and sends the expt initiation message to the eye tracker:
    if n_trials_done == 0:
        # - Update the instruction screen text...
        instr_text_stim.setText(instr_text)
        instr_text_stim.draw()
        flip_time = window.flip()

        # wait until a space key event occurs after the instructions are displayed
        # event.waitKeys(' ')
        event.waitKeys(keyList=['space'])
        # space_key = event.getKeys(keyList=['space'])

    n_trials_done += 1
    print('======TRIAL#' + str(n_trials_done) + '======')

    # Resetting variables:
    this_scue_pitch_hi = None
    scue_onset_frame = None
    scue_delay = None
    this_targ_soa = None
    shutters_shut = None
    shutters_opened = None
    shutter_closing_time = None
    shutter_dur = None

    ## Randomizing variables and assigning the conditions:
    cue_delay = np.random.randint(cue_delay_min, cue_delay_max + 1) / 1000
    if measure:
        cond_str = ''
        this_targ_soa = 600 / 1000
    else:

        # Box (blink) cue delay:
        print('blink cue delay: ' + str(cue_delay))

        # Target location:
        this_targ_loc = trial['targ_right'] * 2 - 1  # converting from binary to [-1, 1]
        if this_targ_loc > 0:
            print('target location: Right')
        else:
            print('target location: Left')
        targ.pos = (targ_off_x * this_targ_loc, 0)

        # Target SOA:
        if exp_name == 'eb1':
            this_targ_soa = 0
        elif exp_name == 'eb2':
            this_targ_soa = trial['soa'] / 1000  # taking SOA from conf-file and converting to seconds (from ms)
        elif exp_name == 'eb3' or exp_name == 'eb4':
            this_targ_soa = soa / 1000
        print('target SOA: ' + str(this_targ_soa) + ' s')

        # Cue validity:
        if trial['cue_valid']:
            print('valid cue')
            if exp_name == 'eb4':
                bcue_box.pos = (0, 0)
                if this_targ_loc > 0:
                    this_scue_pitch_hi = True  # right is higher-pitch
                else:
                    this_scue_pitch_hi = False  # left is lower-pitch
            else:
                bcue_box.pos = (targ_off_x * this_targ_loc, 0)
        else:
            print('invalid cue')
            if exp_name == 'eb4':
                bcue_box.pos = (0, 0)
                if this_targ_loc > 0:
                    this_scue_pitch_hi = False  # left is lower-pitch
                else:
                    this_scue_pitch_hi = True  # right is higher-pitch
            else:
                bcue_box.pos = (-targ_off_x * this_targ_loc, 0)
        print('scue pitch high? = ' + str(this_scue_pitch_hi))  # debug

        # Spatial cue delay (eb4 only):
        if exp_name == 'eb4':
            scue_delay = trial['scue_delay']
            print('spatial cue delay: ' + str(scue_delay))

        noblink_spacebar_recorded = False
        if shutters:
            # noinspection PyUnboundLocalVariable
            shutter_dur = np.random.normal(blink_dur_ave, blink_dur_std)
            print('shutter_dur = ' + str(shutter_dur))
            # if shutter_dur > blink_time_window:
            #     shutter_dur = blink_time_window - 0.11
            shutter_frames = int(shutter_dur * frame_rate)
            shutters_shut = False  # this keeps track of whether the shutters were shut
            shutters_opened = False  # this is an additional tracker for when shutters were already closed this trial
        else:
            shutter_dur = 0

        # Condition string, to pass to the eye tracker, just in case:
        cond_str = ('latency=%s targ_right=%s cue_valid=%s' % (soa, trial['targ_right'], trial['cue_valid']))

    ## Starting the eye-tracking recording.

    # We pump 150 ms delay to allow sufficient time to initiate trials, during which the fixation cross is displayed:
    # flip_time = window.flip()
    fix_cross.draw()

    if not dummy_mode:
        # Starting the recording:
        # take the tracker offline
        tracker.setOfflineMode()
        # tracker.sendMessage('PRE_PUMP %.2f' % flip_time)
        # pylink.pumpDelay(50)
        pylink.msecDelay(50)  # 2020-10-14: replaced the previous line as part of psypy2 to 3 migration

        # Send the standard "TRIALID" message to mark the start of a trial
        # [see Data Viewer User Manual, Section 7: Protocol for EyeLink Data to Viewer Integration]
        tracker.sendMessage('TRIALID %02d' % n_trials_done)

        # Record_status_message : show some info on the host PC
        tracker.sendCommand("record_status_message 'Condition: %s'" % cond_str)

        # Drift check:
        if drift_check:
            # noinspection PyBroadException
            try:
                err = tracker.doDriftCorrect(dr[0] / 2, dr[1] / 2, 1, 1)
            except:
                tracker.doTrackerSetup()
            # read out calibration/drift-correction results:
            print('drift summary = "' + tracker.getCalibrationMessage() + '"')

        # Start recording, parameters specify whether events and samples are stored in file and available over the link
        error = tracker.startRecording(1, 1, 1, 1)
        # pylink.pumpDelay(100)  # wait for 100 ms to make sure data of interest is recorded
        pylink.beginRealTimeMode(100)  # 2020-10-14: replaced the previous line as part of psypy2 to 3 migration

    ## Cycling through the trial phases:
    flip_time = window.flip()
    trial_t_start = flip_time
    if not dummy_mode:
        # Send trial initiation message only post-pump delay:
        if debug:
            print('DEBUG: flip_time is: %.2f' % flip_time)
        tracker.sendMessage('TRIAL_START %.2f' % flip_time)

    # Fixation cross:
    fix_1_frames = int(cue_delay * frame_rate)
    for fix_1_frame in range(fix_1_frames):
        flip_time = frame_routine()

    # ----------------------------------------------
    # Pre-targ frame loop

    bcue_frames = range(int(cue_dur * frame_rate))
    if not measure:
        pretarg_frames = range(int((scue_delay + this_targ_soa) * frame_rate))
        scue_onset_frame = int(scue_delay * frame_rate)
        if debug:
            print('DEBUG: scue onset frame = ' + str(scue_onset_frame))
        ## 2020-10-14: sound prescheduling doesn't seem to work well
        # now = ptb.GetSecs()
        # when_time = now + scue_delay
        # if debug:
        #     print('DEBUG: now time is %.2f' % now)
        #     print('DEBUG: when time is %.2f' % when_time)
        # if this_scue_pitch_hi:
        #     scue_hi.play(when=when_time)
        # else:
        #     scue_lo.play(when=when_time)
        # if debug:
        #     print('DEBUG: Scheduled the sound replay.')
    else:
        pretarg_frames = range(int((1.5 + this_targ_soa) * frame_rate))
    scue_playing = False

    # The location/blink cue:
    event.clearEvents()
    cue_rt_start = flip_time
    cue_rt = 0
    if not dummy_mode:
        tracker.sendMessage('CUE_ONSET %.2f' % flip_time)

    for pretarg_frame in pretarg_frames:

        flip_time = frame_routine()

        if pretarg_frame in bcue_frames:
            bcue_box.draw()

        if not scue_playing and not measure:
            if pretarg_frame >= scue_onset_frame:
                scue_playing = True
                if not dummy_mode:
                    tracker.sendMessage('SCUE_ONSET %.2f' % flip_time)
                if debug:
                    print('DEBUG: spatial cue latency error = ' + str(scue_onset_frame - pretarg_frame))
                    print('DEBUG: audio tone latency = ' + str(flip_time))  # debug
                ## Running without prescheduling:
                if this_scue_pitch_hi:
                    scue_hi.play()
                else:
                    scue_lo.play()

        # In both no-blink and shutters condition, registering the key response to cue:
        if noblink or shutters:
            if not cue_rt:
                cue_rt = monitor_cue_resp(flip_time, cue_rt_start)

        # If the shutters are being used,
        if shutters:
            # ...but are not yet closed...
            if not shutters_shut and not shutters_opened:
                # ...monitoring for the key response to cue...
                if cue_rt > 0:
                    # ...to close the shutters:
                    tracker.sendMessage('CUE_RESP_TIME %.2f' % flip_time)
                    tracker.sendMessage('CUE_RT %.2f' % cue_rt)
                    print('Sent cue RT to the eye tracker.')
                    shutters_shut = close_shutters(dummy_mode)
                    shutter_closing_time = flip_time
            # If the shutters are already closed...
            if shutters_shut and not shutters_opened:
                # ... monitoring the elapsed time from key press to open the shutters when shutter duration is exceeded:
                print(flip_time - shutter_closing_time)
                if (flip_time - shutter_closing_time) >= shutter_dur:
                    shutters_shut = open_shutters(dummy_mode)
                    shutters_opened = True
                    if debug:
                        print(flip_time)
                        print(shutter_closing_time)
                        print(shutter_dur)

        if noblink:
            if not noblink_spacebar_recorded:
                if cue_rt > 0:
                    tracker.sendMessage('CUE_RESP_TIME %.2f' % flip_time)
                    tracker.sendMessage('CUE_RT %.2f' % cue_rt)
                    noblink_spacebar_recorded = True
                    print('Sent cue RT to the eye tracker.')

    # ----------------------------------------------
    # Post-targ frame loop
    if debug:
        print('Post-target frame loop.')
    if not dummy_mode:
        # Passing arbitrary flip times for BLINK LATENCY and WINDOW as they are not relevant
        tracker.sendMessage('BLINK_LATENCY_ONSET %.2f' % flip_time)
        tracker.sendMessage('BLINK_WINDOW_ONSET %.2f' % flip_time)
        tracker.sendMessage('RESPONSE_ONSET %.2f' % flip_time)

        # If the goggles are still somehow shut, opening them forcefully:
        if shutters:
            print('cue_rt ' + str(cue_rt))
            if debug:
                print(shutters_shut)
            if shutters_shut:
                shutters_shut = open_shutters(dummy_mode)
                print('Opened the goggles "forcefully".')
                if debug:
                    print(flip_time)
                    print(shutter_closing_time)
                    print(shutter_dur)

    ## Behavioural response: measuring the reaction time:
    event.clearEvents()

    if not measure:
        if debug:
            print('Tracking response to target.')

        # In the [eb2] experiment, introducing a delay period after which the target is displayed
        if exp_name == 'eb2' and this_targ_soa > 0:

            # # Initiating the delay by flipping once:
            # flip_time = frame_routine()

            # Waiting for the specified time:
            wait(this_targ_soa)

        # Important to measure the flip time for accurate RT measurement just below:
        flip_time = frame_routine()
        if debug:
            print('DEBUG: target latency = ' + str(flip_time))  # debug

        # Trial components pertaining to behavioural response:
        targ_resp_given = False
        rt_start = flip_time

        # Displaying the target and measuring the reaction time.
        while not targ_resp_given:

            # Measuring the time it takes for the behavioural response:
            flip_time = frame_routine()

            # Measuring time elapsed since the start of the trial:
            trial_elapsed_t = flip_time - trial_t_start

            # Drawing the target:
            targ.draw()

            ## Monitoring for key presses:
            arrow_keys = event.getKeys(keyList=['left', 'right'])
            if len(arrow_keys) > 0:
                if 'left' in arrow_keys:
                    print('subject response: Left')
                    beh_resp = -1
                    targ_resp_given = True
                else:
                    print('subject response: Right')
                    beh_resp = 1
                    targ_resp_given = True
                if targ_resp_given:  # this is overwritten every time any key is pressed
                    rt = flip_time - rt_start
                    # noinspection PyUnboundLocalVariable
                    if beh_resp == this_targ_loc:
                        corr_resp = 1  # correct location response
                        accuracy_feedback = 'Correct!'
                    else:
                        corr_resp = 0  # incorrect location response
                        accuracy_feedback = 'INCORRECT!'
                    print('RT=%.2f correct?=%s' % (rt, corr_resp))
                    if not dummy_mode:
                        tracker.sendMessage('TRIAL_RESPONSE %.2f' % flip_time)

    ## Post-trial RT and accuracy
    print('post-trial phase')
    window.flip()
    if not measure:
        # noinspection PyUnboundLocalVariable
        instr_text_stim.setText('Target Location: ' + accuracy_feedback +
                                '\nReaction Time: %.2f' % rt +
                                '\n\nPress the spacebar to continue')
    else:
        instr_text_stim.setText('Press the spacebar to continue')
    instr_text_stim.draw()
    window.flip()

    # wait until a space key event occurs after the instructions are displayed
    event.waitKeys(keyList=['space'])

    flip_time = window.flip()
    if not dummy_mode:
        tracker.sendMessage('TRIAL_END %.2f' % flip_time)

    ## Recording the data
    # noinspection PyUnboundLocalVariable
    if not measure:
        # noinspection PyUnboundLocalVariable
        output_mat[n_trials_done - 1] = {'exp_name': exp_name,
                                         'cue_pred': cue_pred,
                                         'subj': exp_info['subj'],
                                         'cond': exp_info['cond'],
                                         'block': exp_info['block'],
                                         'trial_id': n_trials_done,
                                         'cue_delay': cue_delay,
                                         'targ_right': trial['targ_right'],
                                         'scue_delay': trial['scue_delay'],
                                         'cue_valid': trial['cue_valid'],
                                         'targ_soa': this_targ_soa,
                                         'shutter_dur': shutter_dur,
                                         'trial_start': trial_t_start,
                                         'trial_end': flip_time,
                                         'cue_rt': cue_rt,
                                         'corr_resp': corr_resp,
                                         'rt': rt}
    else:
        output_mat[n_trials_done - 1] = {'exp_name': exp_name,
                                         'subj': exp_info['subj'],
                                         'cond': exp_info['cond'],
                                         'block': exp_info['block'],
                                         'trial_id': n_trials_done,
                                         'cue_delay': cue_delay,
                                         'trial_start': trial_t_start,
                                         'trial_end': flip_time}

    if not dummy_mode:
        ## Stopping the eye tracking
        # send trial variables for Data Viewer integration
        # [see Data Viewer User Manual, Section 7: Protocol for EyeLink Data to Viewer Integration]
        tracker.sendMessage('!V TRIAL_VAR task %s' % cond_str)

        # send a message to mark the end of trial
        # [see Data Viewer User Manual, Section 7: Protocol for EyeLink Data to Viewer Integration]
        tracker.sendMessage('TRIAL_RESULT')
        pylink.endRealTimeMode()
        pylink.pumpDelay(100)
        tracker.stopRecording()

# Finishing the experiment
exit_routine()
