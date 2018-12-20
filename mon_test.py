from psychopy import monitors, visual

dr = (1152, 864)  # display resolution in px
ds = 58  # distance to screen in cm
dd = (34.4, 25.8)  # display dimensions in cm

mon = monitors.Monitor('Station3', width=dd[0], distance=ds)
mon.setSizePix(dr)
window = visual.Window(dr, fullscr=True, monitor=mon, color=[0, 0, 0], units='deg',
                       allowStencil=True, autoLog=False, screen=12, waitBlanking=False)

frame_rate = window.getActualFrameRate()
print('frame rate: ' + str(frame_rate))
