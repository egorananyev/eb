from psychopy import prefs
prefs.hardware['audioLib'] = ['PTB']
import psychtoolbox as ptb
from psychopy import sound, visual, event, monitors

print(sound.Sound)
my_sound = sound.Sound('A')
print('loaded sound')

# Imports associated with the eye tracker:
import pylink
from EyeLinkCoreGraphicsPsychoPy import EyeLinkCoreGraphicsPsychoPy

tracker = pylink.EyeLink('100.1.1.1')

mon = monitors.Monitor('station3')
dr = (1152, 864)
window = visual.Window(dr, fullscr=True, monitor=mon, color='white', units='deg',
                       allowStencil=True, autoLog=False, screen=0, waitBlanking=False)
  
custom_calibration = EyeLinkCoreGraphicsPsychoPy(tracker, window) 
pylink.openGraphicsEx(custom_calibration) 


