import sys, os
from numpy import *

def my_sigint_handler(f,g):
    global frame_forced, frame
    print "Press k to kill, anything else to continue"
    while True:
        s = raw_input("Enter command (k/c/n): ")
        if s=='k':
            print "Exiting"
            sys.exit()
        elif s=='c':
            print "Continuing"
            frame_forced = frame - 2
            if frame_forced<0:
                frame_forced = 0
            break
        elif s=='n':
            frame_forced = int(raw_input("Enter frame: "))
            print "Continuing from frame %d" % frame_forced
            break
        else:
            print "Unknown answer \"%s\"" % s

def existz(afile):
    if os.path.exists(afile):
        return afile
    afilez = afile + ".gz"
    if os.path.exists(afilez):
        return afilez
    return None

## Returns the contenst `a' of `filein', perhaps processed
## with function `fun' as in `a = fun(a,data)' where
## a are the contents of the file. The process is accelerated by
## storing the value in a cached file. 
def read_file_cached(filein,fun=None,data=None,force=False):
    assert os.path.exists(filein),"can't open %s" % filein
    t0 = os.path.getmtime(filein)
    npfile = filein + ".npz"
    np_ok = os.path.exists(npfile)
    if np_ok:
        t1 = os.path.getmtime(npfile)
        np_ok = t1>t0
    if force:
        np_ok = False
    if not np_ok:
        print "reading file %s" % filein
        a = loadtxt(filein)
        print "saving on file %s" % npfile
        if fun is not None:
            a = fun(a,data)
        if isinstance(a,ndarray):
            savez(npfile,a)
        elif isinstance(a,dict):
            savez(npfile,**a)
        else:
            raise TypeError('can''t handle type: ' + type(a))
    print "reading file %s" % npfile
    a = load(npfile)
    if len(a.files)==1 and a.files[0]=='arr_0':
        return a['arr_0']
    else:
        return a

def dump(x):
    savetxt(sys.stdout,x,fmt='%.6e')

def get_cached(key,cache,m,makefun,data=None):
    for q in cache:
        if q[0] == key:
            ## print "key %s, value is cached" % key
            return q[1]
    val = makefun(key,data)
    cache.insert(0,[key,val])
    if len(cache)>m:
        cache.pop()
    return val

def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is one of "yes" or "no".
    """
    valid = {"yes":True,   "y":True,  "ye":True,
             "no":False,     "n":False}
    if default == None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "\
                             "(or 'y' or 'n').\n")

def remove_dir_query(dire):
    if os.path.isdir(dire):
        if query_yes_no("./frames dir exist, delete? ","no"):
            os.system("/bin/rm -rf "+dire)
        else:
            print "ABORT!"
            sys.exit()
    os.mkdir(dire)

def trace(j):
    print "--- TRACE %s ---" % j
