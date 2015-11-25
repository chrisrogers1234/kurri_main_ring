import os
import sys
import shutil
import time
import subprocess
import xboa.Common as Common


PROCESS_LIST = []

def process_wait():
    global PROCESS_LIST
    temp_process_list = []
    for proc in PROCESS_LIST:
        proc.wait()
        print proc.pid, 'returns', proc.returncode
    PROCESS_LIST = []

def process(energy):
    global PROCESS_LIST
    here = os.getcwd()
    dir_name = 'runs/OpalRingTest_ke='+str(energy).rjust(2, '0')
    if os.path.isdir(dir_name):
        shutil.rmtree(dir_name)
    os.makedirs(dir_name)
    time.sleep(1)
    opal_deck = dir_name+'/OpalRingTest.in'
    opal_log = open(dir_name+'/OpalRingTest.log', 'w')
    opal_path = '/home/cr67/OPAL/source/opal_dev/src/trunk/src/opal'
    Common.substitute('OpalRingTest.in',
                      opal_deck,
                      {'__energy__':energy/1000.})
    print 'Running', [opal_path, opal_deck]
    os.chdir(dir_name)
    os.symlink('../../disttest.dat', 'disttest.dat')
    os.symlink('../../fdf-tosca-field-map.table', 'fdf-tosca-field-map.table')
    PROCESS_LIST.append(subprocess.Popen([opal_path, 'OpalRingTest.in'],
                        stdout=opal_log, stderr=subprocess.STDOUT))
    os.chdir(here)

def main():
    process(7.25)
    process(7.5)
    process(7.75)
    process_wait()
    for i in range(14, 17):
        process(i)
    process_wait()
    for i in range(8, 11):
        process(i)
    process_wait()
    for i in range(5, 8):
        process(i)
    process_wait()

if __name__ == "__main__":
    main()
    raw_input()
