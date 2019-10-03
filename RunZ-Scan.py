# This is the main program of DAQ board that incorporates all desired functions. Any further modifications can be built as classes and appended to this.
# noinspection PyUnresolvedReferences
import numpy as np
import time as tm
import daqx as daq
import daqxh as daqh
import sys
# from PyQt5.QtCore import *
# from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
import matplotlib.cm as cm



class Board:    # Hardware control. Connect to the device in the beginning and close the connection only as the last step. This helps avoid initialization problems.

    def __init__(self):
        self.devName = daq.GetDeviceList()[0]        # Find 1st IOtech device

    def open(self): # Connects to the hardware
        print("Connecting to %s\n\n" % self.devName)
        handle = daq.Open(self.devName)          # Open device
        if handle == -1:
            print("Cannot connect to device\n")
            print("Exit")
            sys.exit(handle)
        return handle

    def close(self):
        daq.Close(handle)
        print("\nConnection to %s closed" % self.devName)

    def Ctrig(self):
        # Set port C to output and initialize [C7 ... C0] as LOW
        daq.SetOption(handle,2,daqh.DcofChannel,daqh.DcotP2Local8Mode,daqh.DcovDigitalOutput)
        daq.IOWrite(handle, daqh.DiodtP2Local8, daqh.DiodpP2Local8, 2, daqh.DioepP2, 0)
        # provide hardware trigger by changing C0 from 0 to 1 and back
        # connect TTLTRG to C0
        daq.IOWrite(handle, daqh.DiodtP2Local8, daqh.DiodpP2Local8, 2, daqh.DioepP2, 1)
        daq.IOWrite(handle, daqh.DiodtP2Local8, daqh.DiodpP2Local8, 2, daqh.DioepP2, 0)
        print("Hardware trigger passed.\n")


class Form(QDialog):        # Bring the piezo to focus before performing a z-scan
                            # Contains options for setting by voltage or distance
    def __init__(self, parent=None):

        self.maxVolt = 10.0
        self.minVolt = -10.0
        self.voltage = 1
        self.nscan = 1

        super(Form, self).__init__(parent)

        VoltageLabel = QLabel(" Voltage:         ")
        self.voltageSpinBox = QDoubleSpinBox()
        self.voltageSpinBox.setRange(1.4, 4)
        self.voltageSpinBox.setValue(2)
        self.voltageSpinBox.setSuffix("     ")
        self.voltageSpinBox.setSingleStep(0.01)
        self.voltageSpinBox.setSingleStep(0.01)

        DistanceLabel = QLabel(" Distance:           ")
        self.distanceSpinBox = QDoubleSpinBox()
        self.distanceSpinBox.setRange(15, 45)
        self.distanceSpinBox.setValue(30)
        self.distanceSpinBox.setSuffix("     ")
        self.distanceSpinBox.setSingleStep(1)

        NscanLabel = QLabel(" N Scans:           ")
        self.nscanSpinBox = QDoubleSpinBox()
        self.nscanSpinBox.setRange(1, 30)
        self.nscanSpinBox.setValue(1)
        self.nscanSpinBox.setSuffix("     ")
        self.nscanSpinBox.setSingleStep(1)

        grid = QGridLayout()
        grid.addWidget(VoltageLabel, 0, 0)
        grid.addWidget(self.voltageSpinBox, 0, 1)
        grid.addWidget(DistanceLabel, 1, 0)
        grid.addWidget(self.distanceSpinBox, 1, 1)
        grid.addWidget(NscanLabel, 2, 0)
        grid.addWidget(self.nscanSpinBox, 2, 1)
        self.setLayout(grid)

        self.voltageSpinBox.valueChanged.connect(self.updateVt)
        self.distanceSpinBox.valueChanged.connect(self.updateDt)
        self.nscanSpinBox.valueChanged.connect(self.updateNt)
        # self.connect(self.voltageSpinBox,
        #              SIGNAL("valueChanged(double)"), self.updateVt)
        # self.connect(self.distanceSpinBox,
        #              SIGNAL("valueChanged(double)"), self.updateDt)

        self.setWindowTitle("Piezo focus")
        self.updateDt()

    def updateVt(self):
        self.voltage = self.voltageSpinBox.value()
        self.distanceSpinBox.setValue(self.voltage*15)      # Assumes that piezo moves 15 um per volt.
        cnt_vlt = int(round( (self.voltage-self.minVolt)*65535/(self.maxVolt-self.minVolt) ))
        daq.DacSetOutputMode(handle, daqh.DddtLocal, 0, daqh.DdomVoltage)
        daq.DacWt(handle, daqh.DddtLocal, 0, cnt_vlt)

    def updateDt(self):
        distance = self.distanceSpinBox.value()
        self.voltageSpinBox.setValue(distance/15)           # Assumes that piezo moves 15 um per volt.
        self.voltage = self.voltageSpinBox.value()
        cnt_vlt = int(round( (self.voltage-self.minVolt)*65535/(self.maxVolt-self.minVolt) ))
        daq.DacSetOutputMode(handle, daqh.DddtLocal, 0, daqh.DdomVoltage)
        daq.DacWt(handle, daqh.DddtLocal, 0, cnt_vlt)

    def updateNt(self):
        self.nscan = round(self.nscanSpinBox.value())

class hold(QDialog): #Collect data at a constant voltage

    def __init__(self, parent=None):

        self.holdFREQ=20000
        self.maxVolt = 10.0
        self.minVolt = -10.0
        self.voltage = 2
        self.time=60
        self.holdlength=self.time*self.holdFREQ
        super(hold, self).__init__(parent)

        VoltageLabel = QLabel(" Voltage:         ")
        self.voltageSpinBox = QDoubleSpinBox()
        self.voltageSpinBox.setRange(0.5, 5)
        self.voltageSpinBox.setValue(2)
        self.voltageSpinBox.setSuffix("     ")
        self.voltageSpinBox.setSingleStep(0.01)


        self.startholdbutton = QPushButton("Begin Hold")
        self.startholdbutton.clicked.connect(self.runhold)
        self.stopholdbutton = QPushButton("Stop Hold")
        self.stopholdbutton.clicked.connect(self.stophold)
        TimeLabel = QLabel(" Time:         ")
        self.TimeSpinBox = QDoubleSpinBox()
        self.TimeSpinBox.setRange(1, 300)
        self.TimeSpinBox.setValue(60)
        self.TimeSpinBox.setSuffix("     ")
        self.TimeSpinBox.setSingleStep(1)

        grid = QGridLayout()
        grid.addWidget(VoltageLabel, 0, 0)
        grid.addWidget(TimeLabel, 1, 0)
        grid.addWidget(self.voltageSpinBox, 0, 1)
        grid.addWidget(self.TimeSpinBox, 1, 1)
        grid.addWidget(self.startholdbutton,2,0)
        grid.addWidget(self.stopholdbutton, 3, 0)
        self.setLayout(grid)



        self.setWindowTitle("Hold focus")

        self.voltageSpinBox.valueChanged.connect(self.updateVh)
        self.TimeSpinBox.valueChanged.connect(self.updateTh)

        self.updateVh()
        self.updateTh()


    def updateVh(self): #Update the voltage output from DAQ board
        self.voltage = self.voltageSpinBox.value()
        cnt_vlt = int(round( (self.voltage-self.minVolt)*65535/(self.maxVolt-self.minVolt) ))
        daq.DacSetOutputMode(handle, daqh.DddtLocal, 0, daqh.DdomVoltage)
        daq.DacWt(handle, daqh.DddtLocal, 0, cnt_vlt)
    def updateTh(self): #Change the length of data collection
        self.time = self.TimeSpinBox.value()
        self.holdlength = round(self.time * self.holdFREQ)


    def runhold(self): #Begin data acquisition
        global cellnum
        a1,a2=daq.AdcTransferGetStat(handle)
        if a1 != 0:
            print("Stopping current process")
            daq.AdcDisarm(handle)
            daq.AdcTransferStop(handle)
            daq.DacWaveDisarm(handle, daqh.DddtLocal)


        self.updateVh()
        self.updateTh()
        self.CHANCOUNT = 2
        self.channels = [0,1]                                    # 16 bit counter, 16 bit counter
        self.gains = [daqh.DgainDbd3kX1, daqh.DgainDbd3kX1]  # ignored
        self.flags = [daqh.DafCtr16, daqh.DafCtr16]

        # get read buffer for photon counts
        self.readbuffer    = np.ones(round(self.holdlength*self.CHANCOUNT), dtype=np.uint16)

        # set start and stop conditions
        self.STARTSOURCE	= daqh.DatsExternalTTL
        self.STOPSOURCE	= daqh.DatsScanCount

        daq.AdcDisarm(handle)
        daq.AdcSetAcq(handle, daqh.DaamNShot, 0, self.holdlength)

        # Scan settings
        daq.AdcSetScan(handle, self.channels, self.gains, self.flags)
        # set scan rate
        daq.AdcSetFreq(handle, self.holdFREQ)
        # Setup Channels 0 and 1 (photon counts) for count mode 16 bit, clear on read

        daq.SetOption(handle, self.channels[0], daqh.DcofChannel, daqh.DcotCounterEnhMeasurementMode,
                      daqh.DcovCounterEnhMode_Counter + daqh.DcovCounterEnhCounter_ClearOnRead)
        daq.SetOption(handle, self.channels[1], daqh.DcofChannel, daqh.DcotCounterEnhMeasurementMode,
                      daqh.DcovCounterEnhMode_Counter + daqh.DcovCounterEnhCounter_ClearOnRead)
        # Set buffer location, size and flag settings
        daq.AdcTransferSetBuffer(handle, self.readbuffer, self.holdlength, self.CHANCOUNT,
                                 daqh.DatmUpdateSingle + daqh.DatmCycleOff)
        # Set to Trigger on hardware trigger
        for ch in range(self.CHANCOUNT):
            daq.SetTriggerEvent(handle, self.STARTSOURCE, daqh.DetsRisingEdge, self.channels[ch], self.gains[ch],
                                self.flags[ch], daqh.DaqTypeCounterLocal, 0, 0, daqh.DaqStartEvent)
            # Set to Stop when the requested number of scans is completed
            daq.SetTriggerEvent(handle, self.STOPSOURCE, daqh.DetsRisingEdge, self.channels[ch], self.gains[ch],
                                self.flags[ch], daqh.DaqTypeCounterLocal, 0, 0, daqh.DaqStopEvent)
        daq.AdcTransferStart(handle)
        daq.AdcArm(handle)

        self.win = pg.GraphicsWindow()
        self.win.setWindowTitle('Photon Counts')

        self.p1 = self.win.addPlot()
        self.p1.setLabel('bottom', 'Time', 's')
        self.curve1 = self.p1.plot()
        self.curve2 = self.p1.plot( pen='g')
        self.curvetime = np.arange(self.holdlength) * (self.time / self.holdlength)

        self.timer = pg.QtCore.QTimer()

        self.timer.setInterval(self.time)  # T milliseconds
        self.timer.timeout.connect(self.tichold)

        handle1.Ctrig()
        print("Collecting Data...\n")
        self.timer.start()

    def tichold(self): #Time the data acquisition
        active, retCount = daq.AdcTransferGetStat(handle)
        bufCTR0 = np.array(self.readbuffer.reshape(-1, self.CHANCOUNT)[:, 0])
        bufCTR1 = np.array(self.readbuffer.reshape(-1, self.CHANCOUNT)[:, 1])
        factorval=200
        #print((20*bufCTR0[:(retCount // factorval) * factorval].reshape(-1, factorval).mean(axis=1))[-1])
        self.curve1.setData(self.curvetime[:((retCount // factorval) * factorval):factorval], 20*bufCTR0[:(retCount // factorval) * factorval].reshape(-1, factorval).mean(axis=1))
        self.curve2.setData(self.curvetime[:((retCount // factorval) * factorval):factorval],
                            20 * bufCTR1[:(retCount // factorval) * factorval].reshape(-1, factorval).mean(axis=1))

        #bufCTR0 = self.readbuffer.reshape(-1, self.CHANCOUNT)[:, 0]
        if retCount >= self.holdlength:  # For some reason if retCount is self.SCANS: doesn't work ???
            self.timer.stop()
            # Disarm when completed - placed inside this to avoid possible disarming even before acquiring desired data
            daq.AdcDisarm(handle)
            #daq.DacWaveDisarm(handle, daqh.DddtLocal)
            print("Hold Completed\n")
            # buffer data format [CTR0, CTR1, CTR0, CTR1, CTR0, CTR1, ...]
            bufCTR0 = self.readbuffer.reshape(-1, self.CHANCOUNT)[:, 0]
            bufCTR1 = self.readbuffer.reshape(-1, self.CHANCOUNT)[:, 1]
            #print("bufCTR0 (photon counts)", bufCTR0)
            #print("bufCTR1 (SYNC   counts)", bufCTR1)

            with open('Z:\XTRA folders\Xtra.Jared\Z-Scan DAQ Board\Output\hold cell'+str(cellnum)+'.dat', "wb") as f:
                f.write(bytes(np.column_stack((bufCTR0, bufCTR1))))
            # np.savetxt("PhotonCounts.csv" ,np.column_stack((bufCTR0,bufCTR1,(self.DACwave-self.resolution_half)*self.dV_DAC),delimiter=",",header="Photon Counts,SYNC,Voltage",comments=" ")



    def stophold(self): #Force stop data acquisition
        #try:
            # print('trying')
            # print(type(self))
            # timer=self.timer
            # print('tryingproceed')

        if not hasattr(self, 'timer'):
            print('No hold running\n')
        else:
            self.timer.stop()
            daq.AdcDisarm(handle)
            daq.AdcTransferStop(handle)
            print("Hold Canceled\n")
            active, retCount = daq.AdcTransferGetStat(handle)
            bufCTR0 = self.readbuffer.reshape(-1, self.CHANCOUNT)[:retCount, 0]
            bufCTR1 = self.readbuffer.reshape(-1, self.CHANCOUNT)[:retCount, 1]
            # print("bufCTR0 (photon counts)", bufCTR0)
            # print("bufCTR1 (SYNC   counts)", bufCTR1)
            with open('Z:\XTRA folders\Xtra.Jared\Z-Scan DAQ Board\Output\hold cell'+str(cellnum)+'.dat', "wb") as f:
                f.write(bytes(np.column_stack((bufCTR0,bufCTR1))))

class Zscan:    # Developing the z-scan program with a continuous voltage ramp

    def __init__(self):
        self.i = 0      # Used in self.run() to check if the device is setup at least once

    def setup(self):
        global M2
        # DAC:
        #   Use channel 0 for z-axis output
        #   Use PortAB for digital output streaming (provide a sync pulse indicating the start of the scan)
        #   run DAC output and ADC input (counters) synchronously at a low frequency (traditionally we used a binning frequency of 50000Hz/80 ~ 625Hz
        # dz/dV slope of piezo controller
        if M2 == 1:
            self.dz_dV = 10.0515 # um/V (depends on controller and is an approximate value)
        else:
            self.dz_dV = 15.0773
        self.z0    = 15.0 # um; z-offset (has to be positive to not damage the piezo)
        self.delz  = 24.123602 # um; typical travel of a z-scan
        self.T = 10  # s;  period
        self.N= form.nscan
        self.Thalf = self.T/2
        self.resolution = 1 << 16       # =2**16 (shift bit to left); 16-bit DAC
        self.resolution_half = 1 << 15  # =2**15
        self.Vrange     = 20.0    # V; +-10V output range of DAC
        self.dV_DAC     = self.Vrange/self.resolution
        self.vz    = self.delz/self.Thalf  # speed of the ramp um/s

        # calculate the minimum update rate of DAC to ensure that voltage changes of dV_DAC result in an output update
        self.delV = self.delz/self.dz_dV # voltage amplitude of ramp
        self.V0   = form.voltage   # voltage offset
        self.freq_min = self.delV/self.dV_DAC/self.Thalf # minimum update frequency
        self.FREQ = 250  # set a lower floor for freq of 250 Hz
        # make sure the freq is compatible (commensurable) with the number of data point per scan
        self.nwavehalf = int(round(self.FREQ*self.Thalf)) # number of data points in half of waveform
        self.nwave = self.nwavehalf * 2
        self.FREQ = self.nwavehalf/self.Thalf

        # generate scan waveform
        self.Vtup    = (np.ones(self.nwavehalf).cumsum()-1)/(self.nwavehalf-1)*self.delV     # up-ramp voltage waveform
        self.Vtdown  = np.fliplr([self.Vtup])[0]      # down-ramp voltage waveform
        self.l = round(len(self.Vtup)/2)  # facilitates to make a ramp that starts from middle rather than bottom or top
        self.Vtmid = self.Vtup[self.l]

        self.Vtadjust = self.V0 - self.Vtmid      # facilitates the mid-ramp to start at the focus point set by piezo
        self.Vt = np.concatenate([self.Vtup + self.Vtadjust, self.Vtdown + self.Vtadjust])
        self.Vt = np.tile(self.Vt,self.N)

        # convert Vt into DAC count waveform (-10V = 0,  0V = 32768,  +9.99969482421875V = 65535)
        # ((65535 - 32768 - 65535) + count) * dV_DAC = (count - 32768)*dV_DAC
        self.DACwave = np.around(self.Vt/self.dV_DAC + self.resolution_half).astype('uint16')


        # prepare waveform for digital output; used to mark the start of the scan
        # Note: make sure that the output port is initialized to LOW before start of waveform
        self.COUNT = self.N*self.nwave   # = SCANS
        self.bufSYNC = np.zeros(self.COUNT , dtype=np.uint16)
        self.bufSYNC[0] = 1                                  # |^|________ single pulse at start of digital waveform
        # to read SYNC connect A0 to CNT1

        # prepare ADC readSCAN
        self.SCANS = self.N*self.nwave # = COUNT
        self.CHANCOUNT = 2
        self.channels = [0,1]                                    # 16 bit counter, 16 bit counter
        self.gains    = [daqh.DgainDbd3kX1, daqh.DgainDbd3kX1]   # ignored
        self.flags    = [daqh.DafCtr16,daqh.DafCtr16]

        # get read buffer for photoncounts and sync pulses
        self.readbuffer    = np.ones((self.SCANS*self.CHANCOUNT,), dtype=np.uint16)

        # set start and stop conditions of readSCAN
        self.STARTSOURCE	= daqh.DatsExternalTTL
        self.STOPSOURCE	= daqh.DatsScanCount

        print("Setting up ADC scan...\n")
        daq.AdcSetAcq(handle, daqh.DaamNShot, 0, self.SCANS)
        # Scan settings
        print('m1')
        daq.AdcSetScan(handle, self.channels, self.gains, self.flags)
        # set scan rate
        daq.AdcSetFreq(handle, self.FREQ)
        print('m2')
        # Setup Channel 0 and 1 a(photon counts) for count mode 16 bit, clear on read
        daq.SetOption(handle, self.channels[0], daqh.DcofChannel, daqh.DcotCounterEnhMeasurementMode,
                      daqh.DcovCounterEnhMode_Counter + daqh.DcovCounterEnhCounter_ClearOnRead)
        daq.SetOption(handle, self.channels[1], daqh.DcofChannel, daqh.DcotCounterEnhMeasurementMode,
                      daqh.DcovCounterEnhMode_Counter + daqh.DcovCounterEnhCounter_ClearOnRead)
        print('m3')
        # Set buffer location, size and flag settings
        daq.AdcTransferSetBuffer(handle, self.readbuffer, self.SCANS, self.CHANCOUNT,
                                 daqh.DatmUpdateSingle + daqh.DatmCycleOff)
        print('setup done')
        # Set to Trigger on hardware trigger
        for ch in range(self.CHANCOUNT):
            daq.SetTriggerEvent(handle, self.STARTSOURCE, daqh.DetsRisingEdge, self.channels[ch], self.gains[ch],
                                self.flags[ch], daqh.DaqTypeCounterLocal, 0, 0, daqh.DaqStartEvent)
            # Set to Stop when the requested number of scans is completed
            daq.SetTriggerEvent(handle, self.STOPSOURCE, daqh.DetsRisingEdge, self.channels[ch], self.gains[ch],
                                self.flags[ch], daqh.DaqTypeCounterLocal, 0, 0, daqh.DaqStopEvent)
        print("Setup complete.\n")

    def run(self):
        global cellnum
        a1, a2 = daq.AdcTransferGetStat(handle)
        a1,a2=daq.AdcTransferGetStat(handle)
        if a1 != 0:
            print("Stopping current process")
            daq.AdcDisarm(handle)
            daq.AdcTransferStop(handle)
            daq.DacWaveDisarm(handle, daqh.DddtLocal)
            if hasattr(zscan, 'timer'):
                zscan.stopzscan()
            if hasattr(self, 'timer'):
                self.stophold()

        if self.i is 0:     # Check to see if the zscan is setup at least once.
            print("You need to setup the program first. Starting setup.. \n")
            self.setup()
            #self.i = 1
            print("Starting the process now.. \n")
        # if self.n is not form.nscan:         # Check to see if the number of scans has changed
        #     self.n = form.nscan
        #     self.

        if self.V0 is not form.voltage:         # Check to see if the focus has changed
            self.V0 = form.voltage              # Change the mid-ramp settings caused by change in focus
            self.Vtadjust = self.V0 - self.Vtmid      # facilitates the mid-ramp to start at the focus point set by piezo
            self.Vt = np.concatenate([self.Vtup[self.l:] + self.Vtadjust, self.Vtdown + self.Vtadjust, self.Vtup[:self.l] + self.Vtadjust])
            self.Vt = np.tile(self.Vt,self.N)
            self.DACwave = np.around(self.Vt/self.dV_DAC + self.resolution_half).astype('uint16')

        # initialize DAC output
        # ch0 = V(t)
        daq.DacSetOutputMode     (handle, daqh.DddtLocal, 0, daqh.DdomStaticWave)
        daq.DacWaveSetTrig       (handle, daqh.DddtLocal, 0, daqh.DdtsImmediate, 0)
        daq.DacWaveSetClockSource(handle, daqh.DddtLocal, 0, daqh.DdcsAdcClock)         #slave DAC to ADC clock (synchronous mode)
        daq.DacWaveSetFreq       (handle, daqh.DddtLocal, 0, self.FREQ)                      #set to same frequency as ADC (not sure if this step is necessary)
        daq.DacWaveSetMode       (handle, daqh.DddtLocal, 0, daqh.DdwmInfinite, 0)
        daq.DacWaveSetBuffer     (handle, daqh.DddtLocal, 0, self.DACwave, self.COUNT, daqh.DdtmUserBuffer)
        # ch1 = SYNC(t)
        daq.DacSetOutputMode     (handle, daqh.DddtLocalDigital, 0, daqh.DdomStaticWave)
        daq.DacWaveSetTrig       (handle, daqh.DddtLocalDigital, 0, daqh.DdtsImmediate, 0)
        daq.DacWaveSetClockSource(handle, daqh.DddtLocalDigital, 0, daqh.DdcsAdcClock)
        daq.DacWaveSetFreq       (handle, daqh.DddtLocalDigital, 0, self.FREQ)
        daq.DacWaveSetMode       (handle, daqh.DddtLocalDigital, 0, daqh.DdwmInfinite, 0)
        daq.DacWaveSetBuffer     (handle, daqh.DddtLocalDigital, 0, self.bufSYNC, self.COUNT, daqh.DdtmUserBuffer)



        # begin data acquisition
        daq.AdcTransferStart(handle)
        daq.AdcArm(handle)
        print("waiting for Scan trigger...\n")
        daq.DacWaveArm(handle, daqh.DddtLocal)  #need to arm ADC before DAC to ensure that both wait for hardware trigger

        tm.sleep(1)

        self.win = pg.GraphicsWindow()
        self.win.setWindowTitle('Photon Counts')

        proxy2 = QtGui.QGraphicsProxyWidget()
        stopbutton = QtGui.QPushButton('Stop Entirely')
        proxy2.setWidget(stopbutton)
        stopbutton.clicked.connect(self.stopzscan)

        self.p1 = self.win.addPlot()
        self.p1.setLabel('bottom', 'Time', 's')
        self.curve11 = self.p1.plot()
        self.curve12 = self.p1.plot(pen='g')
        self.win.nextRow()
        self.p2 = self.win.addPlot()

        self.win.nextRow()
        self.p2.setLabel('bottom', 'Voltage', 'V')
        self.curve21 = self.p2.plot()
        self.curve22 = self.p2.plot(pen='g')

        self.win.nextRow()
        self.win.addItem(proxy2)
        self.curvetime = np.arange(self.SCANS)*(self.T*self.N/self.SCANS)       # Used to plot real time data with x-axis as time in seconds

        self.timer = pg.QtCore.QTimer()
        self.timer.setInterval(self.T)      # T milliseconds
        self.timer.timeout.connect(self.update)

        handle1.Ctrig()
        print("Scanning...\n")


        self.timer.start()



    def update(self):
        active, retCount = daq.AdcTransferGetStat(handle)
        bufCTR0 = self.readbuffer.reshape(-1,self.CHANCOUNT)[:,0]
        bufCTR1 = self.readbuffer.reshape(-1, self.CHANCOUNT)[:, 1]
        self.curve11.setData(self.curvetime[:retCount],bufCTR0[:retCount])
        self.curve12.setData(self.curvetime[:retCount], bufCTR1[:retCount])
        self.curve21.setData(self.Vt[:retCount],bufCTR0[:retCount])
        self.curve22.setData(self.Vt[:retCount], bufCTR1[:retCount])
        if retCount >= self.SCANS:          # For some reason if retCount is self.SCANS: doesn't work ???
            self.timer.stop()
            # Disarm when completed - placed inside this to avoid possible disarming even before acquiring desired data
            daq.AdcDisarm(handle)
            daq.DacWaveDisarm(handle, daqh.DddtLocal)
            print("Scan Completed\n")

            cnt_vlt = int(round((self.V0 - form.minVolt) * 65535 / (form.maxVolt - form.minVolt)))
            daq.DacSetOutputMode(handle, daqh.DddtLocal, 0, daqh.DdomVoltage)
            daq.DacWt(handle, daqh.DddtLocal, 0, cnt_vlt)
            # buffer data format [CTR0, CTR1, CTR0, CTR1, CTR0, CTR1, ...]
            bufCTR0 = self.readbuffer.reshape(-1,self.CHANCOUNT)[:,0]
            bufCTR1 = self.readbuffer.reshape(-1,self.CHANCOUNT)[:,1]

            print("bufCTR0 (detector 0 counts)", bufCTR0)
            print("bufCTR1 (detector 1  counts)", bufCTR1)
            np.savetxt("Z:\XTRA folders\Xtra.Jared\Z-Scan DAQ Board\Output\ZScans cell"+str(cellnum)+".csv",np.column_stack((bufCTR0,bufCTR1, ((self.DACwave - self.resolution_half) * self.dV_DAC)*self.dz_dV, (self.DACwave - self.resolution_half) * self.dV_DAC)),delimiter=",")
            # with open('zscan.dat', "wb") as f:
            #     f.write(bytes(np.concatenate([bufCTR0,(self.DACwave - self.resolution_half) * self.dV_DAC))])))
            # np.savetxt("PhotonCounts.csv" ,np.column_stack((bufCTR0,bufCTR1,(self.DACwave-self.resolution_half)*self.dV_DAC),delimiter=",",header="Photon Counts,SYNC,Voltage",comments=" ")
    def stopzscan(self):

        if not hasattr(self, 'timer'):
            print('No scan running\n')
        else:
            self.timer.stop()
            daq.AdcDisarm(handle)
            daq.AdcTransferStop(handle)
            daq.DacWaveDisarm(handle, daqh.DddtLocal)
            print("Scan Canceled\n")

            cnt_vlt = int(round((self.V0 - form.minVolt) * 65535 / (form.maxVolt - form.minVolt)))
            daq.DacSetOutputMode(handle, daqh.DddtLocal, 0, daqh.DdomVoltage)
            daq.DacWt(handle, daqh.DddtLocal, 0, cnt_vlt)
            active, retCount = daq.AdcTransferGetStat(handle)
            bufCTR0 = self.readbuffer.reshape(-1, self.CHANCOUNT)[:retCount, 0]
            bufCTR1 = self.readbuffer.reshape(-1, self.CHANCOUNT)[:retCount, 1]
            # print("bufCTR0 (photon counts)", bufCTR0)
            # print("bufCTR1 (SYNC   counts)", bufCTR1)
            np.savetxt("Z:\XTRA folders\Xtra.Jared\Z-Scan DAQ Board\Output\ZScans cell" + str(cellnum) + ".csv",
                       np.column_stack((bufCTR0,bufCTR1, ((self.DACwave[:retCount] - self.resolution_half) * self.dV_DAC) * self.dz_dV,
                                        (self.DACwave[:retCount] - self.resolution_half) * self.dV_DAC)), delimiter=",")


class ScanThreeHold:    # This is for three-layer geometry brightness deconvolution via alpha matrix formulation
                        # The intensity profile is collected, top and bottom peaks selected, and sequential holds at peaks and center are performed

    def __init__(self):
        self.i = 0      # Used in self.run() to check if the device is setup at least once
        self.numholds=0
    def setup(self):
        global cellnum,M2
        self.numholds = 0
        # DAC:
        #   Use channel 0 for z-axis output
        #   Use PortAB for digital output streaming (provide a sync pulse indicating the start of the scan)
        #   run DAC output and ADC input (counters) synchronously at a low frequency (traditionally we used a binning frequency of 50000Hz/80 ~ 625Hz
        # dz/dV slope of piezo controller
        if M2 == 1:
            self.dz_dV = 10.0515 # um/V (depends on controller and is an approximate value)
        else:
            self.dz_dV = 15.0773
        self.z0    = 15.0 # um; z-offset (has to be positive to not damage the piezo)
        self.delz  = 15 # um; typical travel of a z-scan
        self.T = 5  # s;  period
        self.N= 0.5
        self.Thalf = self.T/2
        self.resolution = 1 << 16       # =2**16 (shift bit to left); 16-bit DAC
        self.resolution_half = 1 << 15  # =2**15
        self.Vrange     = 20.0    # V; +-10V output range of DAC
        self.dV_DAC     = self.Vrange/self.resolution
        self.vz    = self.delz/self.Thalf  # speed of the ramp um/s

        # calculate the minimum update rate of DAC to ensure that voltage changes of dV_DAC result in an output update
        self.delV = self.delz/self.dz_dV # voltage amplitude of ramp
        self.V0   = form.voltage   # voltage offset
        self.freq_min = self.delV/self.dV_DAC/self.Thalf # minimum update frequency
        self.FREQ = 250  # set a lower floor for freq of 1kHz
        # make sure the freq is compatible (commensurable) with the number of data point per scan
        self.nwavehalf = int(round(self.FREQ*self.Thalf)) # number of data points in half of waveform
        self.nwave = self.nwavehalf * 2
        self.FREQ = self.nwavehalf/self.Thalf

        # generate scan waveform
        self.Vtup    = (np.ones(self.nwavehalf).cumsum()-1)/(self.nwavehalf-1)*self.delV     # up-ramp voltage waveform
        self.Vtdown  = np.fliplr([self.Vtup])[0]      # down-ramp voltage waveform
        self.l = round(len(self.Vtup)/2)  # facilitates to make a ramp that starts from middle rather than bottom or top
        self.Vbot = self.Vtup[self.l]
        self.Vtadjust = self.V0 - self.Vbot      # facilitates the mid-ramp to start at the bottom of the scan
        self.Vt = self.Vtup + self.Vtadjust
        # convert Vt into DAC count waveform (-10V = 0,  0V = 32768,  +9.99969482421875V = 65535)
        # ((65535 - 32768 - 65535) + count) * dV_DAC = (count - 32768)*dV_DAC
        self.DACwave = np.around(self.Vt/self.dV_DAC + self.resolution_half).astype('uint16')

        # prepare waveform for digital output; used to mark the start of the scan
        # Note: make sure that the output port is initialized to LOW before start of waveform

        #print(self.COUNT)
        #self.COUNT = self.nwavehalf   # = SCANS
        #print(self.count)
        self.bufSYNC = np.zeros(self.nwavehalf , dtype=np.uint16)
        self.bufSYNC[0] = 1                                  # |^|________ single pulse at start of digital waveform
        # to read SYNC connect A0 to CNT1

        # prepare ADC readSCAN
        #self.SCANS = self.N*self.nwave # = COUNT
        self.CHANCOUNT = 2
        self.channels = [0,1]                                    # 16 bit counter, 16 bit counter
        self.gains    = [daqh.DgainDbd3kX1, daqh.DgainDbd3kX1]   # ignored
        self.flags    = [daqh.DafCtr16,daqh.DafCtr16]

        # get read buffer for photoncounts and sync pulses
        self.readbuffer    = np.ones((self.nwavehalf*self.CHANCOUNT,), dtype=np.uint16)

        # set start and stop conditions of readSCAN
        self.STARTSOURCE	= daqh.DatsExternalTTL
        self.STOPSOURCE	= daqh.DatsScanCount

        #print("Setting up ADC scan...\n")
        daq.AdcSetAcq(handle, daqh.DaamNShot, 0, self.nwavehalf)
        # Scan settings
        daq.AdcSetScan(handle, self.channels, self.gains, self.flags)
        # set scan rate
        daq.AdcSetFreq(handle, self.FREQ)
        # Setup Channel 0 (photon counts) for count mode 16 bit, clear on read
        daq.SetOption(handle, self.channels[0], daqh.DcofChannel, daqh.DcotCounterEnhMeasurementMode,
                      daqh.DcovCounterEnhMode_Counter + daqh.DcovCounterEnhCounter_ClearOnRead)
        daq.SetOption(handle, self.channels[1], daqh.DcofChannel, daqh.DcotCounterEnhMeasurementMode,
                      daqh.DcovCounterEnhMode_Counter + daqh.DcovCounterEnhCounter_ClearOnRead)
        # Set buffer location, size and flag settings
        daq.AdcTransferSetBuffer(handle, self.readbuffer, self.nwavehalf, self.CHANCOUNT,
                                 daqh.DatmUpdateSingle + daqh.DatmCycleOff)
        # Set to Trigger on hardware trigger
        for ch in range(self.CHANCOUNT):
            daq.SetTriggerEvent(handle, self.STARTSOURCE, daqh.DetsRisingEdge, self.channels[ch], self.gains[ch],
                                self.flags[ch], daqh.DaqTypeCounterLocal, 0, 0, daqh.DaqStartEvent)
            # Set to Stop when the requested number of scans is completed
            daq.SetTriggerEvent(handle, self.STOPSOURCE, daqh.DetsRisingEdge, self.channels[ch], self.gains[ch],
                                self.flags[ch], daqh.DaqTypeCounterLocal, 0, 0, daqh.DaqStopEvent)

    def run(self):
        self.setup()



        # if self.V0 is not form.voltage:         # Check to see if the focus has changed
        #     self.V0 = form.voltage              # Change the mid-ramp settings caused by change in focus
        #     self.Vtadjust = self.V0 - self.Vbot      # facilitates the mid-ramp to start at the focus point set by piezo
        #     print(self.Vtadjust)
        #     self.Vt = np.concatenate([self.Vtup[self.l:] + self.Vtadjust, self.Vtdown + self.Vtadjust, self.Vtup[:self.l] + self.Vtadjust])
        #     self.Vt = np.tile(self.Vt,self.N)
        #     self.DACwave = np.around(self.Vt/self.dV_DAC + self.resolution_half).astype('uint16')

        # initialize DAC output
        # ch0 = V(t)
        daq.DacSetOutputMode     (handle, daqh.DddtLocal, 0, daqh.DdomStaticWave)
        daq.DacWaveSetTrig       (handle, daqh.DddtLocal, 0, daqh.DdtsImmediate, 0)
        daq.DacWaveSetClockSource(handle, daqh.DddtLocal, 0, daqh.DdcsAdcClock)         #slave DAC to ADC clock (synchronous mode)
        daq.DacWaveSetFreq       (handle, daqh.DddtLocal, 0, self.FREQ)                      #set to same frequency as ADC (not sure if this step is necessary)
        daq.DacWaveSetMode       (handle, daqh.DddtLocal, 0, daqh.DdwmInfinite, 0)
        daq.DacWaveSetBuffer     (handle, daqh.DddtLocal, 0, self.DACwave, self.nwavehalf, daqh.DdtmUserBuffer)
        # ch1 = SYNC(t)
        daq.DacSetOutputMode     (handle, daqh.DddtLocalDigital, 0, daqh.DdomStaticWave)
        daq.DacWaveSetTrig       (handle, daqh.DddtLocalDigital, 0, daqh.DdtsImmediate, 0)
        daq.DacWaveSetClockSource(handle, daqh.DddtLocalDigital, 0, daqh.DdcsAdcClock)
        daq.DacWaveSetFreq       (handle, daqh.DddtLocalDigital, 0, self.FREQ)
        daq.DacWaveSetMode       (handle, daqh.DddtLocalDigital, 0, daqh.DdwmInfinite, 0)
        daq.DacWaveSetBuffer     (handle, daqh.DddtLocalDigital, 0, self.bufSYNC, self.nwavehalf, daqh.DdtmUserBuffer)



        # begin data acquisition
        daq.AdcTransferStart(handle)
        daq.AdcArm(handle)
        print("waiting for Scan trigger...\n")
        daq.DacWaveArm(handle, daqh.DddtLocal)  #need to arm ADC before DAC to ensure that both wait for hardware trigger

        tm.sleep(1)


        #
        # self.win = pg.GraphicsWindow()
        # self.win.setWindowTitle('Photon Counts')
        #
        # self.p1 = self.win.addPlot()
        # self.p1.setLabel('bottom', 'Time', 's')
        # self.curve1 = self.p1.plot()
        # self.win.nextRow()
        # self.p2 = self.win.addPlot()
        # self.p2.setLabel('bottom', 'Voltage', 'V')
        # self.curve2 = self.p2.plot()
        self.curvetime = np.arange(self.nwavehalf)*(self.T*self.N/self.nwavehalf)       # Used to plot real time data with x-axis as time in seconds

        self.timer = pg.QtCore.QTimer()
        self.timer.setInterval(self.T)      # T milliseconds
        self.timer.timeout.connect(self.update)

        handle1.Ctrig()
        print("Scanning...\n")


        self.timer.start()



    def update(self):
        active, retCount = daq.AdcTransferGetStat(handle)
        # bufCTR0 = self.readbuffer.reshape(-1, self.CHANCOUNT)[:, 0]
        # self.curve1.setData(self.curvetime[:retCount], bufCTR0[:retCount])
        # self.curve2.setData(self.Vt[:retCount], bufCTR0[:retCount])
        if retCount >= self.nwavehalf:          # For some reason if retCount is self.SCANS: doesn't work ???
            self.timer.stop()


            #self.w.scene().sigMouseMoved.connect(mouseClick)

            # Disarm when completed - placed inside this to avoid possible disarming even before acquiring desired data
            daq.AdcDisarm(handle)
            daq.DacWaveDisarm(handle, daqh.DddtLocal)


            # buffer data format [CTR0, CTR1, CTR0, CTR1, CTR0, CTR1, ...]
            bufCTR0 = self.readbuffer.reshape(-1,self.CHANCOUNT)[:,0]
            bufCTR1 = self.readbuffer.reshape(-1,self.CHANCOUNT)[:,1]
            self.bufCTR0=bufCTR0
            #print("bufCTR0 (photon counts)", bufCTR0)
            #print("bufCTR1 (SYNC   counts)", bufCTR1)
            #print(self.cellnum)
            np.savetxt("Z:\XTRA folders\Xtra.Jared\Z-Scan DAQ Board\Output\Fast Z-Scan cell"+str(cellnum)+".csv",np.column_stack((bufCTR0, bufCTR1,((self.DACwave - self.resolution_half) * self.dV_DAC)*self.dz_dV, (self.DACwave - self.resolution_half) * self.dV_DAC)),delimiter=",")
            cnt_vlt = int(round((self.Vtup[0] - form.minVolt) * 65535 / (form.maxVolt - form.minVolt)))
            daq.DacSetOutputMode(handle, daqh.DddtLocal, 0, daqh.DdomVoltage)
            daq.DacWt(handle, daqh.DddtLocal, 0, cnt_vlt)
            self.determineholdvoltages()
            # np.savetxt("PhotonCounts.csv" ,np.column_stack((bufCTR0,bufCTR1,(self.DACwave-self.resolution_half)*self.dV_DAC),delimiter=",",header="Photon Counts,SYNC,Voltage",comments=" ")
            #("Scan Completed. Select three locations for holds.\n")

            # def onClick(event):
            #     items = self.win.scene().items(event.scenePos())
            #     print("Plots:", [x for x in items if isinstance(x, pg.PlotItem)])
            #
            # self.win.scene().sigMouseMoved.connect(onClick)

             #print(self.win.mapFromScene(event.scenePos()))

    def determineholdvoltages(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(self.Vt, self.bufCTR0)
        ax.plot(self.Vt, self.bufCTR1,'g')
        self.holdlocs= [0]
        self.holdfs=[]
        def resetlocs(event):
            print('Resetting Locations')
            self.holdlocs=[]

        def onclick(event):
            global ix, iy
            ix, iy = event.xdata, event.ydata
            # print('x = %d, y = %d' % (ix, iy))
            if len(self.holdlocs) >= 1:
                print(ix)

            if ix != None:
                self.holdlocs.append((ix))
                self.holdfs.append((iy))
            if len(self.holdlocs) == 4:
                cnt_vlt = int(round((self.Vt[0] - self.minVolt) * 65535 / (self.maxVolt - self.minVolt)))
                daq.DacSetOutputMode(handle, daqh.DddtLocal, 0, daqh.DdomVoltage)
                daq.DacWt(handle, daqh.DddtLocal, 0, cnt_vlt)
                self.holdlocs=self.holdlocs[1:]
                fig.canvas.mpl_disconnect(cid)
                np.savetxt("Z:\XTRA folders\Xtra.Jared\Z-Scan DAQ Board\Output\Hold Voltages cell"+str(cellnum)+".csv",
                           self.holdlocs,
                           delimiter=",", comments=" ")
                self.runthreeholds()


        axreset = plt.axes([0.15, 0.95, 0.25, 0.05])
        breset = Button(axreset, 'Reset Locations')
        breset.on_clicked(resetlocs)

        cid = fig.canvas.mpl_connect('button_release_event', onclick)

        plt.show()

    def runthreeholds(self):
        self.minVolt=-10
        self.maxVolt=10
        if self.numholds > 0:
            filename = 'Z:\XTRA folders\Xtra.Jared\Z-Scan DAQ Board\Output\hold'+str(self.numholds)+'cell'+str(cellnum)+'.dat'
            with open(filename, "wb") as f:
                f.write(bytes(np.column_stack((self.bufCTR0, self.bufCTR1))))
        # if self.numholds == 2:
        #     filename = 'Z:\XTRA folders\Xtra.Jared\Z-Scan DAQ Board\Output\Second hold cell'+str(cellnum)+'.dat'
        #     with open(filename, "wb") as f:
        #         f.write(bytes(self.bufCTR0))
        # if self.numholds == 3:
        #     filename = 'Z:\XTRA folders\Xtra.Jared\Z-Scan DAQ Board\Output\Third hold cell'+str(cellnum)+'.dat'
        #     with open(filename, "wb") as f:
        #         f.write(bytes(self.bufCTR0))

        if self.numholds < 3:
            self.voltage = self.holdlocs[self.numholds]
            cnt_vlt = int(round((self.voltage - self.minVolt) * 65535 / (self.maxVolt - self.minVolt)))
            daq.DacSetOutputMode(handle, daqh.DddtLocal, 0, daqh.DdomVoltage)
            daq.DacWt(handle, daqh.DddtLocal, 0, cnt_vlt)
            print('Voltage set to',self.voltage)
            print('Running hold\n')
            self.runsinglehold()

        else:
            print('All holds complete, acquiring z-scans')
            zscan.run()


    def runsinglehold(self):
        self.holdFREQ = 20000
        self.time = 60.0
        self.holdlength = round(self.time * self.holdFREQ)
        self.CHANCOUNT = 2
        self.channels = [0, 1]  # 16 bit counter, 16 bit counter
        self.gains = [daqh.DgainDbd3kX1, daqh.DgainDbd3kX1]  # ignored
        self.flags = [daqh.DafCtr16, daqh.DafCtr16]

        # get read buffer for photoncounts and sync pulses
        self.readbuffer = np.ones(round(self.holdlength * self.CHANCOUNT), dtype=np.uint16)

        # set start and stop conditions of readSCAN
        self.STARTSOURCE = daqh.DatsExternalTTL
        self.STOPSOURCE = daqh.DatsScanCount

        daq.AdcDisarm(handle)
        daq.AdcSetAcq(handle, daqh.DaamNShot, 0, self.holdlength)

        # Scan settings
        daq.AdcSetScan(handle, self.channels, self.gains, self.flags)
        # set scan rate
        daq.AdcSetFreq(handle, self.holdFREQ)
        # Setup Channel 0 and 1 (photon counts) for count mode 16 bit, clear on read
        daq.SetOption(handle, self.channels[0], daqh.DcofChannel, daqh.DcotCounterEnhMeasurementMode,
                      daqh.DcovCounterEnhMode_Counter + daqh.DcovCounterEnhCounter_ClearOnRead)
        daq.SetOption(handle, self.channels[1], daqh.DcofChannel, daqh.DcotCounterEnhMeasurementMode,
                      daqh.DcovCounterEnhMode_Counter + daqh.DcovCounterEnhCounter_ClearOnRead)
        # Set buffer location, size and flag settings
        daq.AdcTransferSetBuffer(handle, self.readbuffer, self.holdlength, self.CHANCOUNT,
                                 daqh.DatmUpdateSingle + daqh.DatmCycleOff)
        # Set to Trigger on hardware trigger
        for ch in range(self.CHANCOUNT):
            daq.SetTriggerEvent(handle, self.STARTSOURCE, daqh.DetsRisingEdge, self.channels[ch], self.gains[ch],
                                self.flags[ch], daqh.DaqTypeCounterLocal, 0, 0, daqh.DaqStartEvent)
            # Set to Stop when the requested number of scans is completed
            daq.SetTriggerEvent(handle, self.STOPSOURCE, daqh.DetsRisingEdge, self.channels[ch], self.gains[ch],
                                self.flags[ch], daqh.DaqTypeCounterLocal, 0, 0, daqh.DaqStopEvent)
        daq.AdcTransferStart(handle)
        daq.AdcArm(handle)

        self.win = pg.GraphicsWindow()
        self.win.setWindowTitle('Photon Counts')

        proxy = QtGui.QGraphicsProxyWidget()
        button = QtGui.QPushButton('Restart Hold')
        proxy.setWidget(button)
        button.clicked.connect(self.restarthold)

        proxy2 = QtGui.QGraphicsProxyWidget()
        stopbutton = QtGui.QPushButton('Stop Entirely')
        proxy2.setWidget(stopbutton)
        stopbutton.clicked.connect(self.stopentirely)

        self.p1 = self.win.addPlot()
        self.p1.setLabel('bottom', 'Time', 's')
        self.curve1 = self.p1.plot()
        self.curve2 = self.p1.plot(pen='g')
        self.curvetime = np.arange(self.holdlength) * (self.time / self.holdlength)
        self.win.nextRow()
        self.win.addItem(proxy)
        self.win.addItem(proxy2)
        self.timer = pg.QtCore.QTimer()

        self.timer.setInterval(self.time)  # T milliseconds
        self.timer.timeout.connect(self.tichold)

        handle1.Ctrig()
        print("Collecting Data...\n")
        self.timer.start()

    def tichold(self):
        active, retCount = daq.AdcTransferGetStat(handle)
        bufCTR0 = np.array(self.readbuffer.reshape(-1, self.CHANCOUNT)[:, 0])
        bufCTR1 = np.array(self.readbuffer.reshape(-1, self.CHANCOUNT)[:, 1])
        factorval = 200
        # print((20*bufCTR0[:(retCount // factorval) * factorval].reshape(-1, factorval).mean(axis=1))[-1])
        self.curve1.setData(self.curvetime[:((retCount // factorval) * factorval):factorval],
                            20 * bufCTR0[:(retCount // factorval) * factorval].reshape(-1, factorval).mean(axis=1))
        self.curve2.setData(self.curvetime[:((retCount // factorval) * factorval):factorval],
                            20 * bufCTR1[:(retCount // factorval) * factorval].reshape(-1, factorval).mean(axis=1))

        # bufCTR0 = self.readbuffer.reshape(-1, self.CHANCOUNT)[:, 0]
        if retCount >= self.holdlength:  # For some reason if retCount is self.SCANS: doesn't work ???
            self.timer.stop()
            # Disarm when completed - placed inside this to avoid possible disarming even before acquiring desired data
            daq.AdcDisarm(handle)
            # daq.DacWaveDisarm(handle, daqh.DddtLocal)
            print("Hold Completed\n")
            # buffer data format [CTR0, CTR1, CTR0, CTR1, CTR0, CTR1, ...]
            self.bufCTR0 = self.readbuffer.reshape(-1, self.CHANCOUNT)[:, 0]
            self.bufCTR1 = self.readbuffer.reshape(-1, self.CHANCOUNT)[:, 1]
            self.numholds=self.numholds+1
            # print("bufCTR0 (photon counts)", bufCTR0)
            # print("bufCTR1 (SYNC   counts)", bufCTR1)
            self.runthreeholds()
            # np.savetxt("PhotonCounts.csv" ,np.column_stack((bufCTR0,bufCTR1,(self.DACwave-self.resolution_half)*self.dV_DAC),delimiter=",",header="Photon Counts,SYNC,Voltage",comments=" ")

    def restarthold(self):
        if not hasattr(self, 'timer'):
            print('No hold running\n')
        else:
            self.timer.stop()
            daq.AdcDisarm(handle)
            daq.AdcTransferStop(handle)
            self.runthreeholds()

    def stopentirely(self):
        if hasattr(self, 'timer'):
            print('Stopped\n')
            self.timer.stop()
            daq.AdcDisarm(handle)
            daq.AdcTransferStop(handle)
            daq.DacWaveDisarm(handle, daqh.DddtLocal)

class ScanplusNHolds(QDialog):    # Repeat of above ScanThreeHold but with arbitrary number of holds

    def __init__(self, parent=None):

        self.i = 0      # Used in self.run() to check if the device is setup at least once
        self.numholds=0
        self.holdFREQ=20000
        self.maxVolt = 10.0
        self.minVolt = -10.0
        self.time=60
        self.holdlength=self.time*self.holdFREQ
        super(ScanplusNHolds, self).__init__(parent)

        Nholdslabel = QLabel(" Number of Holds:         ")
        self.NholdsSpinBox = QDoubleSpinBox()
        self.NholdsSpinBox.setRange(2, 10)
        self.NholdsSpinBox.setValue(5)
        self.NholdsSpinBox.setSuffix("     ")
        self.NholdsSpinBox.setSingleStep(1)


        self.prelimscan = QPushButton("Preliminary Scan")
        self.prelimscan.clicked.connect(self.run)

        TimeLabel = QLabel(" Time Per Hold:         ")
        self.TimeSpinBox = QDoubleSpinBox()
        self.TimeSpinBox.setRange(1, 300)
        self.TimeSpinBox.setValue(30)
        self.TimeSpinBox.setSuffix("     ")
        self.TimeSpinBox.setSingleStep(1)

        grid = QGridLayout()
        grid.addWidget(Nholdslabel, 0, 0)
        grid.addWidget(TimeLabel, 1, 0)
        grid.addWidget(self.NholdsSpinBox, 0, 1)
        grid.addWidget(self.TimeSpinBox, 1, 1)
        grid.addWidget(self.prelimscan,2,0)
        self.setLayout(grid)



        self.setWindowTitle("N-Hold Setup")

        self.NholdsSpinBox.valueChanged.connect(self.updateNh)
        self.TimeSpinBox.valueChanged.connect(self.updateTh)

        self.updateNh()
        self.updateTh()


    def updateNh(self):
        self.Nholds = round(self.NholdsSpinBox.value())
    def updateTh(self):
        self.time = self.TimeSpinBox.value()
        self.holdlength = round(self.time * self.holdFREQ)

    def setup(self):
        global cellnum,M2
        self.numholds = 0
        # DAC:
        #   Use channel 0 for z-axis output
        #   Use PortAB for digital output streaming (provide a sync pulse indicating the start of the scan)
        #   run DAC output and ADC input (counters) synchronously at a low frequency (traditionally we used a binning frequency of 50000Hz/80 ~ 625Hz
        # dz/dV slope of piezo controller
        if M2 == 1:
            self.dz_dV = 10.0515 # um/V (depends on controller and is an approximate value)
        else:
            self.dz_dV = 15.0773
        self.z0    = 15.0 # um; z-offset (has to be positive to not damage the piezo)
        self.delz  = 15 # um; typical travel of a z-scan
        self.T = 5  # s;  period
        self.N= 0.5
        self.Thalf = self.T/2
        self.resolution = 1 << 16       # =2**16 (shift bit to left); 16-bit DAC
        self.resolution_half = 1 << 15  # =2**15
        self.Vrange     = 20.0    # V; +-10V output range of DAC
        self.dV_DAC     = self.Vrange/self.resolution
        self.vz    = self.delz/self.Thalf  # speed of the ramp um/s

        # calculate the minimum update rate of DAC to ensure that voltage changes of dV_DAC result in an output update
        self.delV = self.delz/self.dz_dV # voltage amplitude of ramp
        self.V0   = form.voltage   # voltage offset
        self.freq_min = self.delV/self.dV_DAC/self.Thalf # minimum update frequency
        self.FREQ = 250  # set a lower floor for freq of 1kHz
        # make sure the freq is compatible (commensurable) with the number of data point per scan
        self.nwavehalf = int(round(self.FREQ*self.Thalf)) # number of data points in half of waveform
        self.nwave = self.nwavehalf * 2
        self.FREQ = self.nwavehalf/self.Thalf

        # generate scan waveform
        self.Vtup    = (np.ones(self.nwavehalf).cumsum()-1)/(self.nwavehalf-1)*self.delV     # up-ramp voltage waveform
        self.Vtdown  = np.fliplr([self.Vtup])[0]      # down-ramp voltage waveform
        self.l = round(len(self.Vtup)/2)  # facilitates to make a ramp that starts from middle rather than bottom or top
        self.Vbot = self.Vtup[self.l]
        self.Vtadjust = self.V0 - self.Vbot      # facilitates the mid-ramp to start at the bottom of the scan
        self.Vt = self.Vtup + self.Vtadjust
        # convert Vt into DAC count waveform (-10V = 0,  0V = 32768,  +9.99969482421875V = 65535)
        # ((65535 - 32768 - 65535) + count) * dV_DAC = (count - 32768)*dV_DAC
        self.DACwave = np.around(self.Vt/self.dV_DAC + self.resolution_half).astype('uint16')

        # prepare waveform for digital output; used to mark the start of the scan
        # Note: make sure that the output port is initialized to LOW before start of waveform

        #print(self.COUNT)
        #self.COUNT = self.nwavehalf   # = SCANS
        #print(self.count)
        self.bufSYNC = np.zeros(self.nwavehalf , dtype=np.uint16)
        self.bufSYNC[0] = 1                                  # |^|________ single pulse at start of digital waveform
        # to read SYNC connect A0 to CNT1

        # prepare ADC readSCAN
        #self.SCANS = self.N*self.nwave # = COUNT
        self.CHANCOUNT = 2
        self.channels = [0,1]                                    # 16 bit counter, 16 bit counter
        self.gains    = [daqh.DgainDbd3kX1, daqh.DgainDbd3kX1]   # ignored
        self.flags    = [daqh.DafCtr16,daqh.DafCtr16]

        # get read buffer for photoncounts and sync pulses
        self.readbuffer    = np.ones((self.nwavehalf*self.CHANCOUNT,), dtype=np.uint16)

        # set start and stop conditions of readSCAN
        self.STARTSOURCE	= daqh.DatsExternalTTL
        self.STOPSOURCE	= daqh.DatsScanCount

        #print("Setting up ADC scan...\n")
        daq.AdcSetAcq(handle, daqh.DaamNShot, 0, self.nwavehalf)
        # Scan settings
        daq.AdcSetScan(handle, self.channels, self.gains, self.flags)
        # set scan rate
        daq.AdcSetFreq(handle, self.FREQ)
        # Setup Channel 0 and 1 (photon counts) for count mode 16 bit, clear on read
        daq.SetOption(handle, self.channels[0], daqh.DcofChannel, daqh.DcotCounterEnhMeasurementMode,
                      daqh.DcovCounterEnhMode_Counter + daqh.DcovCounterEnhCounter_ClearOnRead)
        daq.SetOption(handle, self.channels[1], daqh.DcofChannel, daqh.DcotCounterEnhMeasurementMode,
                      daqh.DcovCounterEnhMode_Counter + daqh.DcovCounterEnhCounter_ClearOnRead)
        # Set buffer location, size and flag settings
        daq.AdcTransferSetBuffer(handle, self.readbuffer, self.nwavehalf, self.CHANCOUNT,
                                 daqh.DatmUpdateSingle + daqh.DatmCycleOff)
        # Set to Trigger on hardware trigger
        for ch in range(self.CHANCOUNT):
            daq.SetTriggerEvent(handle, self.STARTSOURCE, daqh.DetsRisingEdge, self.channels[ch], self.gains[ch],
                                self.flags[ch], daqh.DaqTypeCounterLocal, 0, 0, daqh.DaqStartEvent)
            # Set to Stop when the requested number of scans is completed
            daq.SetTriggerEvent(handle, self.STOPSOURCE, daqh.DetsRisingEdge, self.channels[ch], self.gains[ch],
                                self.flags[ch], daqh.DaqTypeCounterLocal, 0, 0, daqh.DaqStopEvent)

    def run(self):
        self.setup()


        # if self.V0 is not form.voltage:         # Check to see if the focus has changed
        #     self.V0 = form.voltage              # Change the mid-ramp settings caused by change in focus
        #     self.Vtadjust = self.V0 - self.Vbot      # facilitates the mid-ramp to start at the focus point set by piezo
        #     print(self.Vtadjust)
        #     self.Vt = np.concatenate([self.Vtup[self.l:] + self.Vtadjust, self.Vtdown + self.Vtadjust, self.Vtup[:self.l] + self.Vtadjust])
        #     self.Vt = np.tile(self.Vt,self.N)
        #     self.DACwave = np.around(self.Vt/self.dV_DAC + self.resolution_half).astype('uint16')

        # initialize DAC output
        # ch0 = V(t)
        daq.DacSetOutputMode     (handle, daqh.DddtLocal, 0, daqh.DdomStaticWave)
        daq.DacWaveSetTrig       (handle, daqh.DddtLocal, 0, daqh.DdtsImmediate, 0)
        daq.DacWaveSetClockSource(handle, daqh.DddtLocal, 0, daqh.DdcsAdcClock)         #slave DAC to ADC clock (synchronous mode)
        daq.DacWaveSetFreq       (handle, daqh.DddtLocal, 0, self.FREQ)                      #set to same frequency as ADC (not sure if this step is necessary)
        daq.DacWaveSetMode       (handle, daqh.DddtLocal, 0, daqh.DdwmInfinite, 0)
        daq.DacWaveSetBuffer     (handle, daqh.DddtLocal, 0, self.DACwave, self.nwavehalf, daqh.DdtmUserBuffer)
        # ch1 = SYNC(t)
        daq.DacSetOutputMode     (handle, daqh.DddtLocalDigital, 0, daqh.DdomStaticWave)
        daq.DacWaveSetTrig       (handle, daqh.DddtLocalDigital, 0, daqh.DdtsImmediate, 0)
        daq.DacWaveSetClockSource(handle, daqh.DddtLocalDigital, 0, daqh.DdcsAdcClock)
        daq.DacWaveSetFreq       (handle, daqh.DddtLocalDigital, 0, self.FREQ)
        daq.DacWaveSetMode       (handle, daqh.DddtLocalDigital, 0, daqh.DdwmInfinite, 0)
        daq.DacWaveSetBuffer     (handle, daqh.DddtLocalDigital, 0, self.bufSYNC, self.nwavehalf, daqh.DdtmUserBuffer)



        # begin data acquisition
        daq.AdcTransferStart(handle)
        daq.AdcArm(handle)
        print("waiting for Scan trigger...\n")
        daq.DacWaveArm(handle, daqh.DddtLocal)  #need to arm ADC before DAC to ensure that both wait for hardware trigger

        tm.sleep(1)


        #
        # self.win = pg.GraphicsWindow()
        # self.win.setWindowTitle('Photon Counts')
        #
        # self.p1 = self.win.addPlot()
        # self.p1.setLabel('bottom', 'Time', 's')
        # self.curve1 = self.p1.plot()
        # self.win.nextRow()
        # self.p2 = self.win.addPlot()
        # self.p2.setLabel('bottom', 'Voltage', 'V')
        # self.curve2 = self.p2.plot()
        self.curvetime = np.arange(self.nwavehalf)*(self.T*self.N/self.nwavehalf)       # Used to plot real time data with x-axis as time in seconds

        self.timer = pg.QtCore.QTimer()
        self.timer.setInterval(self.T)      # T milliseconds
        self.timer.timeout.connect(self.updatescan)

        handle1.Ctrig()
        print("Scanning...\n")


        self.timer.start()



    def updatescan(self):
        active, retCount = daq.AdcTransferGetStat(handle)
        # bufCTR0 = self.readbuffer.reshape(-1, self.CHANCOUNT)[:, 0]
        # self.curve1.setData(self.curvetime[:retCount], bufCTR0[:retCount])
        # self.curve2.setData(self.Vt[:retCount], bufCTR0[:retCount])
        if retCount >= self.nwavehalf:          # For some reason if retCount is self.SCANS: doesn't work ???
            self.timer.stop()


            #self.w.scene().sigMouseMoved.connect(mouseClick)

            # Disarm when completed - placed inside this to avoid possible disarming even before acquiring desired data
            daq.AdcDisarm(handle)
            daq.DacWaveDisarm(handle, daqh.DddtLocal)


            # buffer data format [CTR0, CTR1, CTR0, CTR1, CTR0, CTR1, ...]
            bufCTR0 = self.readbuffer.reshape(-1,self.CHANCOUNT)[:,0]
            bufCTR1 = self.readbuffer.reshape(-1,self.CHANCOUNT)[:,1]
            self.bufCTR0=bufCTR0
            self.bufCTR1=bufCTR1
            #print("bufCTR0 (photon counts)", bufCTR0)
            #print("bufCTR1 (SYNC   counts)", bufCTR1)
            #print(self.cellnum)
            np.savetxt("Z:\XTRA folders\Xtra.Jared\Z-Scan DAQ Board\Output\Fast Z-Scan cell"+str(cellnum)+".csv",np.column_stack((bufCTR0, bufCTR1,((self.DACwave - self.resolution_half) * self.dV_DAC)*self.dz_dV, (self.DACwave - self.resolution_half) * self.dV_DAC)),delimiter=",")
            cnt_vlt = int(round((self.Vtup[0] - form.minVolt) * 65535 / (form.maxVolt - form.minVolt)))
            daq.DacSetOutputMode(handle, daqh.DddtLocal, 0, daqh.DdomVoltage)
            daq.DacWt(handle, daqh.DddtLocal, 0, cnt_vlt)
            self.determineholdvoltages()
            # np.savetxt("PhotonCounts.csv" ,np.column_stack((bufCTR0,bufCTR1,(self.DACwave-self.resolution_half)*self.dV_DAC),delimiter=",",header="Photon Counts,SYNC,Voltage",comments=" ")
            #("Scan Completed. Select three locations for holds.\n")

            # def onClick(event):
            #     items = self.win.scene().items(event.scenePos())
            #     print("Plots:", [x for x in items if isinstance(x, pg.PlotItem)])
            #
            # self.win.scene().sigMouseMoved.connect(onClick)

             #print(self.win.mapFromScene(event.scenePos()))

    def determineholdvoltages(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(self.Vt, self.bufCTR0,'k')
        ax.plot(self.Vt, self.bufCTR1,'g')
        self.holdlocs= [0]
        self.holdfs=[]
        def resetlocs(event):
            print('Resetting Locations')
            self.holdlocs=[]

        def onclick(event):
            global ix, iy
            ix, iy = event.xdata, event.ydata
            # print('x = %d, y = %d' % (ix, iy))
            if len(self.holdlocs) >= 1:
                print(ix)

            if ix != None:
                self.holdlocs.append((ix))
                self.holdfs.append((iy))
            if len(self.holdlocs) == 3:
                cnt_vlt = int(round((self.Vt[0] - self.minVolt) * 65535 / (self.maxVolt - self.minVolt)))
                daq.DacSetOutputMode(handle, daqh.DddtLocal, 0, daqh.DdomVoltage)
                daq.DacWt(handle, daqh.DddtLocal, 0, cnt_vlt)
                Nholds=self.Nholds
                self.holdlocs=np.cumsum((self.holdlocs[2]-self.holdlocs[1])/(Nholds-1)*np.ones(Nholds))+self.holdlocs[1]-(self.holdlocs[2]-self.holdlocs[1])/(Nholds-1)
                fig.canvas.mpl_disconnect(cid)
                np.savetxt("Z:\XTRA folders\Xtra.Jared\Z-Scan DAQ Board\Output\Hold Voltages cell"+str(cellnum)+".csv",
                           self.holdlocs,
                           delimiter=",", comments=" ")

                plt.close()
                self.runmultipleholds()




        axreset = plt.axes([0.15, 0.95, 0.25, 0.05])
        breset = Button(axreset, 'Reset Locations')
        breset.on_clicked(resetlocs)

        cid = fig.canvas.mpl_connect('button_release_event', onclick)

        plt.show()

    def runmultipleholds(self):
        self.minVolt=-10
        self.maxVolt=10
        if self.numholds > 0:
            filename = 'Z:\XTRA folders\Xtra.Jared\Z-Scan DAQ Board\Output\hold'+str(self.numholds)+' cell'+str(cellnum)+'.dat'
            with open(filename, "wb") as f:
                f.write(bytes(np.column_stack((self.bufCTR0, self.bufCTR1))))
        if self.numholds < (self.Nholds):
            self.voltage = self.holdlocs[self.numholds]
            cnt_vlt = int(round((self.voltage - self.minVolt) * 65535 / (self.maxVolt - self.minVolt)))
            daq.DacSetOutputMode(handle, daqh.DddtLocal, 0, daqh.DdomVoltage)
            daq.DacWt(handle, daqh.DddtLocal, 0, cnt_vlt)
            print('Voltage set to',self.voltage)
            print('Running hold\n')
            self.runsinglehold()

        else:
            print('All holds complete, acquiring z-scans')
            zscan.run()


    def runsinglehold(self):
        self.CHANCOUNT = 2
        self.channels = [0, 1]  # 16 bit counter, 16 bit counter
        self.gains = [daqh.DgainDbd3kX1, daqh.DgainDbd3kX1]  # ignored
        self.flags = [daqh.DafCtr16, daqh.DafCtr16]

        # get read buffer for photoncounts and sync pulses
        self.readbuffer = np.ones(round(self.holdlength * self.CHANCOUNT), dtype=np.uint16)

        # set start and stop conditions of readSCAN
        self.STARTSOURCE = daqh.DatsExternalTTL
        self.STOPSOURCE = daqh.DatsScanCount

        daq.AdcDisarm(handle)
        daq.AdcSetAcq(handle, daqh.DaamNShot, 0, self.holdlength)

        # Scan settings
        daq.AdcSetScan(handle, self.channels, self.gains, self.flags)
        # set scan rate
        daq.AdcSetFreq(handle, self.holdFREQ)
        # Setup Channel 0 (photon counts) for count mode 16 bit, clear on read
        daq.SetOption(handle, self.channels[0], daqh.DcofChannel, daqh.DcotCounterEnhMeasurementMode,
                      daqh.DcovCounterEnhMode_Counter + daqh.DcovCounterEnhCounter_ClearOnRead)
        daq.SetOption(handle, self.channels[1], daqh.DcofChannel, daqh.DcotCounterEnhMeasurementMode,
                      daqh.DcovCounterEnhMode_Counter + daqh.DcovCounterEnhCounter_ClearOnRead)
        # Set buffer location, size and flag settings
        daq.AdcTransferSetBuffer(handle, self.readbuffer, self.holdlength, self.CHANCOUNT,
                                 daqh.DatmUpdateSingle + daqh.DatmCycleOff)
        # Set to Trigger on hardware trigger
        for ch in range(self.CHANCOUNT):
            daq.SetTriggerEvent(handle, self.STARTSOURCE, daqh.DetsRisingEdge, self.channels[ch], self.gains[ch],
                                self.flags[ch], daqh.DaqTypeCounterLocal, 0, 0, daqh.DaqStartEvent)
            # Set to Stop when the requested number of scans is completed
            daq.SetTriggerEvent(handle, self.STOPSOURCE, daqh.DetsRisingEdge, self.channels[ch], self.gains[ch],
                                self.flags[ch], daqh.DaqTypeCounterLocal, 0, 0, daqh.DaqStopEvent)
        daq.AdcTransferStart(handle)
        daq.AdcArm(handle)

        self.win = pg.GraphicsWindow()
        self.win.setWindowTitle('Photon Counts')

        proxy = QtGui.QGraphicsProxyWidget()
        button = QtGui.QPushButton('Restart Hold')
        proxy.setWidget(button)
        button.clicked.connect(self.restarthold)

        proxy2 = QtGui.QGraphicsProxyWidget()
        stopbutton = QtGui.QPushButton('Stop Entirely')
        proxy2.setWidget(stopbutton)
        stopbutton.clicked.connect(self.stopentirely)

        self.p1 = self.win.addPlot()
        self.p1.setLabel('bottom', 'Time', 's')
        self.curve1 = self.p1.plot()
        self.curve2 = self.p1.plot(pen='g')
        self.curvetime = np.arange(self.holdlength) * (self.time / self.holdlength)
        self.win.nextRow()
        self.win.addItem(proxy)
        self.win.addItem(proxy2)
        self.timer = pg.QtCore.QTimer()

        self.timer.setInterval(self.time)  # T milliseconds
        self.timer.timeout.connect(self.tichold)

        handle1.Ctrig()
        print("Collecting Data...\n")
        self.timer.start()

    def tichold(self):
        active, retCount = daq.AdcTransferGetStat(handle)
        bufCTR0 = np.array(self.readbuffer.reshape(-1, self.CHANCOUNT)[:, 0])
        bufCTR1 = np.array(self.readbuffer.reshape(-1, self.CHANCOUNT)[:, 1])
        factorval = 200
        # print((20*bufCTR0[:(retCount // factorval) * factorval].reshape(-1, factorval).mean(axis=1))[-1])
        self.curve1.setData(self.curvetime[:((retCount // factorval) * factorval):factorval],
                            20 * bufCTR0[:(retCount // factorval) * factorval].reshape(-1, factorval).mean(axis=1))
        self.curve2.setData(self.curvetime[:((retCount // factorval) * factorval):factorval],
                            20 * bufCTR1[:(retCount // factorval) * factorval].reshape(-1, factorval).mean(axis=1))

        # bufCTR0 = self.readbuffer.reshape(-1, self.CHANCOUNT)[:, 0]
        if retCount >= self.holdlength:  # For some reason if retCount is self.SCANS: doesn't work ???
            self.timer.stop()
            # Disarm when completed - placed inside this to avoid possible disarming even before acquiring desired data
            daq.AdcDisarm(handle)
            # daq.DacWaveDisarm(handle, daqh.DddtLocal)
            print("Hold Completed\n")
            # buffer data format [CTR0, CTR1, CTR0, CTR1, CTR0, CTR1, ...]
            self.bufCTR0 = self.readbuffer.reshape(-1, self.CHANCOUNT)[:, 0]
            self.bufCTR1 = self.readbuffer.reshape(-1, self.CHANCOUNT)[:, 1]
            self.numholds=self.numholds+1
            # print("bufCTR0 (photon counts)", bufCTR0)
            # print("bufCTR1 (SYNC   counts)", bufCTR1)
            self.runmultipleholds()
            # np.savetxt("PhotonCounts.csv" ,np.column_stack((bufCTR0,bufCTR1,(self.DACwave-self.resolution_half)*self.dV_DAC),delimiter=",",header="Photon Counts,SYNC,Voltage",comments=" ")

    def restarthold(self):
        if not hasattr(self, 'timer'):
            print('No hold running\n')
        else:
            self.timer.stop()
            daq.AdcDisarm(handle)
            daq.AdcTransferStop(handle)
            self.runmultipleholds()

    def stopentirely(self):
        if hasattr(self, 'timer'):
            print('Stopped\n')
            self.timer.stop()
            daq.AdcDisarm(handle)
            daq.AdcTransferStop(handle)

class MainGUI(QtGui.QDialog):

    def __init__(self, parent = None):

        super(MainGUI, self).__init__(parent)
        self.setWindowTitle("DAQ board")
        cellnumlabel = QLabel("Cell Number:           ")
        self.cellnum=1
        self.cellnumSpinBox = QDoubleSpinBox()
        self.cellnumSpinBox.setRange(1, 1000)
        self.cellnumSpinBox.setValue(1)
        self.cellnumSpinBox.setSuffix("     ")
        self.cellnumSpinBox.setSingleStep(1)
        self.button1 = QPushButton("Setup")
        self.button2 = QPushButton("Run")
        self.button3 = QPushButton("Scan Params")
        self.button4 = QPushButton("Hold")
        self.button5 = QPushButton("Three Point Collection")
        self.Nholdbutton = QPushButton("N Point Collection")
        self.m1check = QCheckBox ('M1')
        self.button1.clicked.connect(self.one)
        self.button2.clicked.connect(self.two)
        self.button3.clicked.connect(self.focus)
        self.button4.clicked.connect(self.hold)
        self.button5.clicked.connect(self.scanhold)
        self.Nholdbutton.clicked.connect(self.Nhold)
        #ZscanLabel = QLabel("Z-scan: ")
        grid = QGridLayout()
        #grid.addWidget(ZscanLabel, 0, 0)
        #grid.addWidget(self.button1, 0, 2)
        grid.addWidget(cellnumlabel,3,0)
        grid.addWidget(self.cellnumSpinBox,3,1)
        grid.addWidget(self.button2, 0, 1)
        grid.addWidget(self.button3, 0, 0)
        grid.addWidget(self.button4, 1, 0)
        grid.addWidget(self.button5, 1, 1)
        grid.addWidget(self.Nholdbutton,2,0)
        grid.addWidget(self.m1check,4,0)
        self.setLayout(grid)
        self.m1check.stateChanged.connect(self.setM1)
        self.cellnumSpinBox.valueChanged.connect(self.updatecellnum)
        self.setM1()
        self.updatecellnum()
    def updatecellnum(self):
        self.cellnum = round(self.cellnumSpinBox.value())
        global cellnum
        cellnum=self.cellnum
    def one(self):
        zscan.setup()
    def two(self):
        zscan.run()
    def focus(self):
        form.show()
    def galvo(self):
        galvo.show()
    def hold(self):
        hold.show()
    def scanhold(self):
        scanthreehold.run()
    def Nhold(self):
        scanplusnholds.show()
    def setM1(self):
        global M2
        if self.m1check.isChecked() == 1:
            # global M2
            M2=0
        else:
            # global M2
            M2=1







if __name__ == "__main__":
    handle1 = Board()
    handle = handle1.open()

    app = QApplication([])

    scanthreehold= ScanThreeHold()
    zscan = Zscan()
    hold=hold()
    form = Form()
    scanplusnholds=ScanplusNHolds()
    GUI = MainGUI()
    GUI.show()
    app.exec_()
    handle1.close()
