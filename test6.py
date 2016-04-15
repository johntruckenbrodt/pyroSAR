
import Queue
from threading import Thread
from subprocess import Popen, PIPE


def enqueue_output(out, queue):
    for line in iter(out.readline, b''):
        queue.put(line)
    out.close()


def getOutput(outQueue):
    outStr = ''
    try:
        while True:  # Adds output from the Queue until it is empty
            outStr += outQueue.get_nowait()
    except Queue.Empty:
        return outStr


p = Popen(['create_dem_par', '/geonfs01_vol1/ve39vem/S1/test_camarque/test.par'], stdin=PIPE, stdout=PIPE, stderr=PIPE, shell=False, universal_newlines=True)

outQueue = Queue.Queue()
errQueue = Queue.Queue()

outThread = Thread(target=enqueue_output, args=(p.stdout, outQueue))
errThread = Thread(target=enqueue_output, args=(p.stderr, errQueue))

outThread.daemon = True
errThread.daemon = True

outThread.start()
errThread.start()

# someInput = raw_input("Input: ")

p.stdin.write("UTM\n")
p.stdin.write("WGS84\n")
p.stdin.write("1\n")
p.stdin.write("32\n")
p.stdin.write("10000000.\n")
p.stdin.write("testdem\n")
p.stdin.write("REAL*4\n")
p.stdin.write("0.0\n")
p.stdin.write("1.0\n")
p.stdin.write("2000\n")
p.stdin.write("4000\n")
p.stdin.write("-20 20\n")
p.stdin.write("432654 34765\n")
errors = getOutput(errQueue)
output = getOutput(outQueue)

#####################


import sys
from subprocess import PIPE, Popen, STDOUT
from threading import Thread

from Queue import Queue, Empty


ON_POSIX = 'posix' in sys.builtin_module_names

def enqueue_output(out, queue):
    for line in iter(out.readline, b''):
        queue.put(line)
    out.close()

p = Popen(['create_dem_par', '/geonfs01_vol1/ve39vem/S1/test_camarque/test.par'], stdin=PIPE, stdout=PIPE, stderr=STDOUT, bufsize=1, close_fds=ON_POSIX)
# p = Popen(['create_dem_par', '/geonfs01_vol1/ve39vem/S1/test_camarque/test.par'], stdin=PIPE, stdout=PIPE, stderr=STDOUT, shell=False, universal_newlines=True)
q = Queue()
t = Thread(target=enqueue_output, args=(p.stdout, q))
t.daemon = True # thread dies with the program
t.start()

# r = Queue()
# s = Thread(target=enqueue_output, args=(p.stderr, r))
# s.daemon = True
# s.start()

q.get_nowait()
# r.get_nowait()
p.stdin.write("UTM\n")