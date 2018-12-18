####!/usr/bin/arch -i386 /usr/bin/python # -*- coding: utf-8 -*-
from __future__ import print_function

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
AllOff = 'x
"""

from __future__ import division  # so that 1/3=0.333 instead of 1/3=0
from psychopy import visual, core, event, gui, sound, monitors
import numpy as np
import os
from psychopy.data import TrialHandler, importConditions
import pandas as pd
from datetime import datetime

# Imports associated with the eye tracker:
import pylink
from EyeLinkCoreGraphicsPsychoPy import EyeLinkCoreGraphicsPsychoPy

## Initial variables.
# experiment modes:
toshi = True
dummy_mode = True
# experiment variables:
exp_name = 'eb1'
trial_n = 15  # trials per row; 15 gives 300 trials
fix_1_dur = .4  # the time frame before the cue, if any
blink_latency_min = 90  # these are in ms, because we need a random _integer_ in this range
blink_latency_max = 350  # Note! the range is actually 240-500 ms, but 150 are already included in eyelink waiting
# the time window for the blink - quite conservative - should include the whole blink, but is independent of
# the blink start/end
blink_time_window = .3
# display dimensions:
if toshi:
    dr = (576, 432)
    ds = 61
    dd = (17.2, 12.9)  # display dimensions in cm
else:
    dr = (1152, 864)  # display resolution in px
    ds = 58  # distance to screen in cm
    dd = (34.4, 25.8)  # display dimensions in cm
# cue:
cue_off_y = .8
cue_dur = .2
beep_dur = .05
# target:
targ_off_x = 3
targ_diam = .5

## getting user info about the experiment session:
exp_info = {u'expt': exp_name, u'subj': u'', u'cond': u'd', u'sess': u''}
# conditions: 't'=training, 'c'=control, 'a'=artificial blink, 'v'=voluntary blink, 'd'=debug
exp_name = exp_info['expt']
dlg = gui.DlgFromDict(dictionary=exp_info, title=exp_name)  # dialogue box
if not dlg.OK:
    core.quit()  # user pressed cancel
exp_info['time'] = datetime.now().strftime('%Y-%m-%d_%H%M')

# Assigning conditions:
debug = False
eye_tracking = True  # true by default
training = False
voluntary = False
shutters = False
print('Condition: ' + exp_info['cond'])
if exp_info['cond'] == 'd':
    debug = True
    eye_tracking = False
if exp_info['cond'] == 't':
    trial_n = 1
    training = True
if exp_info['cond'] == 'v':
    voluntary = True
if exp_info['cond'] == 'a':
    shutters = True
    import serial

    ser = serial.Serial('/dev/ttyACM0', 9600)
    ser.write('c')

## Input and output

# condition file:
exp_conditions = importConditions('cond-files/cond_' + exp_name + '_' + exp_info['cond'] + '.xlsx')
trials = TrialHandler(exp_conditions, trial_n, extraInfo=exp_info)

# output file:
out_file_name = '%s_subj-%s_cond-%s_sess-%s_%s' % (exp_name, exp_info['subj'], exp_info['cond'],
                                                   exp_info['sess'], exp_info['time'])
out_file_path = '..' + os.sep + 'data' + os.sep + out_file_name

# output matrix:
output_mat = {}

## EyeLink setup
if not dummy_mode:
    tracker = pylink.EyeLink('100.1.1.1')
else:
    tracker = pylink.EyeLink(None)

# Note that the file name cannot exceeds 8 characters. Open eyelink data files as early to record as possible.
edf_subj_data_dir_name = out_file_name
edf_data_file_path = '..' + os.sep + 'edf_data' + os.sep + edf_subj_data_dir_name
if not os.path.exists(edf_data_file_path):
    os.makedirs(edf_data_file_path)
edf_data_file_name = datetime.now().strftime('%m%d%H%M') + '.edf'  # to avoid overwriting, naming MMDDHHmm
tracker.openDataFile(edf_data_file_name)
# add personalized data file header (preamble text)
tracker.sendCommand("add_file_preamble_text 'Study: Influence of blinks on attentional cueing'")

## Monitor setup
if toshi:
    mon = monitors.Monitor('Dell', width=dd[0], distance=ds)
    window = visual.Window(size=dr, monitor=mon, fullscr=False, screen=1, units='deg')
else:
    # you MUST specify the physical properties of your monitor first, otherwise you won't be able to properly use
    # different screen "units" in psychopy. One may define his/her monitor object within the GUI, but
    # I find it is a better practice to put things all under control in the experimental script instead.
    mon = monitors.Monitor('Station3', width=dd[0], distance=ds)
    mon.setSizePix(dr)
    window = visual.Window(dr, fullscr=True, monitor=mon, color=[0, 0, 0], units='pix',
                           allowStencil=True, autoLog=False)

if debug:
    frame_rate = 60
else:
    frame_rate = window.getActualFrameRate()
    print('frame rate: ' + str(frame_rate))
    if frame_rate < 100:
        print('WARNING! The measured frame rate is lower than expected')

## Initialize the stimuli and instructions
instruction_text = "Please press the ''space'' key\nto start the experiment"
instructions_text_stim = visual.TextStim(window, text=instruction_text, height=.4)
# todo: measure cross size
fix_cross = visual.TextStim(window, text='+', bold='True', pos=[0, 0], rgb=1, height=.3)

# auditory tone (for blink condition):
print('using %s (with %s) for sounds' % (sound.audioLib, sound.audioDriver))
beep = sound.Sound('A', octave=5, sampleRate=44100, secs=beep_dur, stereo=True)
beep.setVolume(.1)  # this isn't needed in Aaron's paradigm -- the sound is softer

# cue:
arrow_vert = [(.5, 0), (0, .3), (0, .1), (-.5, .1), (-.5, -.1), (0, -.1), (0, -.3)]
cue_arrow = visual.ShapeStim(window, vertices=arrow_vert, fillColor='black', size=.3, lineColor='black',
                             pos=(0, cue_off_y))

# target:
targ = visual.Circle(window, radius=targ_diam / 2, edges=32, pos=(10, 0), fillColor='white')

## Eye-tracking calibration:
# call the custom calibration routine "EyeLinkCoreGraphicsPsychopy.py", instead of the default
# routines that were implemented in SDL
custom_calibration = EyeLinkCoreGraphicsPsychoPy(tracker, window)
pylink.openGraphicsEx(custom_calibration)

## STEP V: Set up the tracker
# we need to put the tracker in offline mode before we change its configrations
tracker.setOfflineMode()
# sampling rate, 250, 500, 1000, or 2000; this command won't work for EyeLInk II/I
tracker.sendCommand('sample_rate 500')
# inform the tracker the resolution of the subject display
# [see Eyelink Installation Guide, Section 8.4: Customizing Your PHYSICAL.INI Settings ]
tracker.sendCommand("screen_pixel_coords = 0 0 %d %d" % (scnWidth - 1, scnHeight - 1))
# save display resolution in EDF data file for Data Viewer integration purposes
# [see Data Viewer User Manual, Section 7: Protocol for EyeLink Data to Viewer Integration]
tracker.sendMessage("DISPLAY_COORDS = 0 0 %d %d" % (scnWidth - 1, scnHeight - 1))
# specify the calibration type, H3, HV3, HV5, HV13 (HV = horizontal/vertical),
tracker.sendCommand("calibration_type = HV9")  # tracker.setCalibrationType('HV9') also works, see the Pylink manual
# the model of the tracker, 1-EyeLink I, 2-EyeLink II, 3-Newer models (100/1000Plus/DUO)
eyelinkVer = tracker.getTrackerVersion()
# turn off 'scene link' camera stuff (EyeLink II/I only)
if eyelinkVer == 2:
    tracker.sendCommand("scene_camera_gazemap = NO")
# Set the tracker to parse Events using "GAZE" (or "HREF") data
tracker.sendCommand("recording_parse_type = GAZE")
# Online parser configuration: 0-> standard/cognitive, 1-> sensitive/psychophysiological
# the Parser for EyeLink I is more conservative, see below
# [see Eyelink User Manual, Section 4.3: EyeLink Parser Configuration]
if eyelinkVer >= 2:
    tracker.sendCommand('select_parser_configuration 0')
# get Host tracking software version
hostVer = 0
if eyelinkVer == 3:
    tvstr = tracker.getTrackerVersionString()
    vindex = tvstr.find("EYELINK CL")
    hostVer = int(float(tvstr[(vindex + len("EYELINK CL")):].strip()))

## Sending eye-tracker commands and calibration:
# specify the EVENT and SAMPLE data that are stored in EDF or retrievable from the Link
# See Section 4 Data Files of the EyeLink user manual
tracker.sendCommand("file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT")
tracker.sendCommand("link_event_filter = LEFT,RIGHT,FIXATION,FIXUPDATE,SACCADE,BLINK,BUTTON,INPUT")
if hostVer >= 4:
    tracker.sendCommand("file_sample_data  = LEFT,RIGHT,GAZE,AREA,GAZERES,STATUS,HTARGET,INPUT")
    tracker.sendCommand("link_sample_data  = LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS,HTARGET,INPUT")
else:
    tracker.sendCommand("file_sample_data  = LEFT,RIGHT,GAZE,AREA,GAZERES,STATUS,INPUT")
    tracker.sendCommand("link_sample_data  = LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS,INPUT")

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


# Also no variation across frames, but only available upon call, which is made only in key registering phase.
def exit_routine():
    if shutters:
        ser.write('z')

    # Behavioural data output:
    data_columns = ['exp_name', 'subj', 'cond', 'sess', 'trial_id', 'targ_right', 'cue_valid',
                    'blink_latency', 'trial_start', 'trial_end', 'corr_resp', 'rt']
    pd.DataFrame.from_dict(output_mat, orient='index').to_csv(out_file_path + '.csv', index=False,
                                                              columns=data_columns)
    # .to_csv(out_file_path + '.csv', index=False, columns=data_columns)
    print('output file path is ' + out_file_path)

    # EDF output:
    # close the EDF data file
    tracker.setOfflineMode()
    tracker.closeDataFile()
    pylink.pumpDelay(50)

    # Get the EDF data and say goodbye
    instructions_text_stim.setText('Recording data...')
    instructions_text_stim.draw()
    window.flip()
    tracker.receiveDataFile(edf_data_file_name, edf_data_file_path + os.sep + edf_data_file_name)

    # close the link to the tracker
    tracker.close()

    # close the graphics
    pylink.closeGraphics()
    window.close()
    core.quit()


## Initiating the trial loop

n_trials_done = 0

for trial in trials:

    ## First trial initiates instructions and sends the expt initiation message to the eye tracker:
    if n_trials_done == 0:
        # - Update the instruction screen text...
        instructions_text_stim.setText(instruction_text)
        instructions_text_stim.draw()
        flip_time = window.flip()

        # wait until a space key event occurs after the instructions are displayed
        event.waitKeys(' ')

    n_trials_done += 1
    print('======TRIAL#' + str(n_trials_done) + '======')

    if debug:
        if n_trials_done > 1:
            # noinspection PyUnboundLocalVariable
            iti_dur = flip_time - iti_end_trial
            print('inter-trial duration: %.3f' % iti_dur)

    ## Randomizing variables and assigning the conditions:
    # Randomize the duration of the post-cue fixation & converting to sec:
    blink_latency = np.random.randint(blink_latency_min,
                                      blink_latency_max + 1) / 1000  # max value has to be one up
    if debug:
        print('blink_latency = %.3f' % blink_latency)

    # Target location:
    # this_targ_loc = np.random.randint(2) * 2 - 1
    this_targ_loc = trial['targ_right'] * 2 - 1
    if this_targ_loc > 0:
        print('target location: Right')
    else:
        print('target location: Left')
    targ.pos = (targ_off_x * this_targ_loc, 0)

    # Cue validity:
    cue_dir = (trial['cue_valid'] * 2 - 1) * this_targ_loc
    # Logic: First, the cue validity is converted from binary [0, 1] to [-1, 1].
    # It is then multiplied by the target location, which is either left [-1] or right [1].
    # E.g., if the cue is valid for a target that appears on the right, cue direction is 1*1, rightward.
    # If the cue is invalid for such a target, cue direction is -1*1, leftward.
    cue_arrow.ori = 90 * (cue_dir - 1)
    # Logic: If cue_dir == 1, 90 * 0 = 0, rightward orientation. If cue_dir == -1, 90 * (-2) = -180, leftward.
    if trial['cue_valid']:
        print('valid cue')
    else:
        print('invalid cue')

    # Condition string, to pass to the eye tracker, just in case:
    cond_str = ('latency=%s targ_right=%s cue_valid=%s' % (blink_latency, trial['targ_right'], trial['cue_valid']))

    ## Starting the eye-tracking recording.

    # We pump 150 ms delay to allow sufficient time to initiate trials, during which the fixation cross is displayed:
    flip_time = window.flip()
    fix_cross.draw()

    # Starting the recording:
    # take the tracker offline
    tracker.setOfflineMode()
    pylink.pumpDelay(50)

    # send the standard "TRIALID" message to mark the start of a trial
    # [see Data Viewer User Manual, Section 7: Protocol for EyeLink Data to Viewer Integration]
    tracker.sendMessage('TRIALID')

    # record_status_message : show some info on the host PC
    tracker.sendCommand("record_status_message 'Condition: %s'" % cond_str)

    # drift check
    try:
        err = tracker.doDriftCorrect(scnWidth / 2, scnHeight / 2, 1, 1)
    except:
        tracker.doTrackerSetup()

    # read out calibration/drift-correction results:
    print(tracker.getCalibrationMessage())

    # start recording, parameters specify whether events and samples are stored in file and available over the link
    error = tracker.startRecording(1, 1, 1, 1)
    pylink.pumpDelay(100)  # wait for 100 ms to make sure data of interest is recorded

    ## Cycling through the trial phases:
    trial_t_start = flip_time
    if debug:
        trial_elapsed_frames = 0  # counting frames for frame skip test

    # Fixation cross:
    fix_1_frames = int(fix_1_dur * frame_rate)
    for fix_1_frame in range(fix_1_frames):
        frame_routine()
        if debug:
            # noinspection PyUnboundLocalVariable
            trial_elapsed_frames += 1

    # The rest of the period without the beep, but with the cue:
    cue_frames = int(cue_dur * frame_rate)
    for cue_frame in range(cue_frames):
        frame_routine()
        cue_arrow.draw()
        if debug:
            # noinspection PyUnboundLocalVariable
            trial_elapsed_frames += 1

    # The brief period with the beep:
    beep_frames = int(beep_dur * frame_rate)
    for beep_frame in range(beep_frames):
        frame_routine()
        beep.play()
        if debug:
            # noinspection PyUnboundLocalVariable
            trial_elapsed_frames += 1

    # Blink latency = the fixation period after the beep:
    blink_latency_frames = int(blink_latency * frame_rate)
    for blink_latency_frame in range(blink_latency_frames):
        flip_time = frame_routine()
        if debug:
            # noinspection PyUnboundLocalVariable
            trial_elapsed_frames += 1

    # Real or simulated blink follow the same timeline:
    blink_time_period_frames = int(blink_time_window * frame_rate)
    if shutters:
        ser.write('z')
        print('Closed the goggles.')
    for blink_time_period_frame in range(blink_time_period_frames):
        flip_time = frame_routine()
        if debug:
            # noinspection PyUnboundLocalVariable
            trial_elapsed_frames += 1
    if shutters:
        ser.write('c')
        print('Opened the goggles.')

    ## Behavioural response: measuring the reaction time:

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

        if debug:
            # noinspection PyUnboundLocalVariable
            trial_elapsed_frames += 1

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
                if beh_resp == this_targ_loc:
                    corr_resp = 1  # correct location response
                else:
                    corr_resp = 0  # incorrect location response
                print('RT=%.2f correct?=%s' % (rt, corr_resp))
                if debug:  # in debug mode, check if the frame rate looks okay
                    # noinspection PyUnboundLocalVariable
                    frame_skip_check(trial_elapsed_t, trial_elapsed_frames)
                    iti_end_trial = flip_time

    ## Recording the data
    # noinspection PyUnboundLocalVariable
    output_mat[n_trials_done - 1] = {'exp_name': exp_name,
                                     'subj': exp_info['subj'],
                                     'cond': exp_info['cond'],
                                     'sess': exp_info['sess'],
                                     'trial_id': n_trials_done,
                                     'targ_right': trial['targ_right'],
                                     'cue_valid': trial['cue_valid'],
                                     'blink_latency': blink_latency,
                                     'trial_start': trial_t_start,
                                     'trial_end': flip_time,
                                     'corr_resp': corr_resp,
                                     'rt': rt}

    ## Stopping the eye tracking
    # send trial variables for Data Viewer integration
    # [see Data Viewer User Manual, Section 7: Protocol for EyeLink Data to Viewer Integration]
    tracker.sendMessage('!V TRIAL_VAR task %s' % cond_str)

    # send a message to mark the end of trial
    # [see Data Viewer User Manual, Section 7: Protocol for EyeLink Data to Viewer Integration]
    tracker.sendMessage('TRIAL_RESULT')
    pylink.pumpDelay(100)
    tracker.stopRecording()

# Finishing the experiment
exit_routine()
