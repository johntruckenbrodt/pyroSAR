##############################################################
# GUI Interface for creating GAMMA process interface windows
# John Truckenbrodt 2014-2015
##############################################################

"""
this interface is used for creating child windows for defining parameters and excuting chosen commands
the arguments needed are passed form the main script pythonland.py in nested lists
by modifying the variable args.extra additional text can be defined for specific command child windows, which will be printed above the parameter selection form
the size of the window is defined according to the parameters needed for the chosen command
if a corresponding parameter file exists, the parameter values of this file are entered into the form, otherwise the default values (passed from class Main of script pythonland.py)
are entered
a button is created for starting the process
"""

import os
import Tkinter
from fileQuery import FileQuery
from ancillary import ReadPar
from auxiliary import Environment, makeform, execute


class Dialog(Tkinter.Toplevel):
    def __init__(self, args):
        Tkinter.Toplevel.__init__(self)
        self.config(Environment.bcolor)
        self.resizable(width=Tkinter.FALSE, height=Tkinter.FALSE)
        Tkinter.Frame(self, Environment.bcolor, height=10, bd=5).pack()

        self.title = args[0][0]

        self.action = args[0] if len(args[0]) == 1 else [args[0][1], args[0][2]]

        self.objList = []

        # define names of dropdown menu options
        self.dropnames = {}

        # extra vertical window space
        self.space = 0

        # add window label from args_extra
        if args[0][0] in Environment.args_extra:
            self.label = Tkinter.Label(self, Environment.header_ops, text=Environment.args_extra[args[0][0]])
            self.label.pack()
            self.space = 10

        if len(args[0]) > 1:
            parname = os.path.join(Environment.workdir.get(), "PAR", os.path.splitext(args[0][2])[0]+".par")
            if os.path.isfile(parname):
                self.params =ReadPar(parname, type="exe")

        for i in range(0, len(args[1])):
            if args[1][i] in ["import directory", "SRTM archive"]:
                self.objList.append(FileQuery(self, args[1][i], 2))
            elif args[1][i] == "output file":
                self.objList.append(FileQuery(self, args[1][i], 3))
            else:
                self.objList.append(FileQuery(self, args[1][i], 1))

        # define y-dimension of the dialog window
        self.ydim = 80 + self.space + len(args[1]) * 30 + len(args[3]) * 30

        # check whether parameter textfile for the chosen command exists;
        # if so then use the parameters from the file, otherwise use the defaults defined in class Main (pythonland.py)

        if hasattr(self, "params"):
            if "SRTM_archive" in self.params.index:
                self.params.index.remove("SRTM_archive")
            self.ents = makeform(self, self.params.index, [getattr(self.params, attr) for attr in self.params.index if attr != "SRTM_archive"])
        else:
            self.ents = makeform(self, args[-2], args[-1])

        # set window appearance
        self.geometry("600x" + str(self.ydim))
        Tkinter.Frame(self, bg="white", height=2).pack({"fill": "x"})
        Tkinter.Frame(self, Environment.bcolor, height=10, bd=5).pack()

        # create execution button
        self.Action = Tkinter.Button(self, Environment.button_ops, text="Action", padx=40, command=lambda: execute(self.action, self.objList, self.ents))
        self.Action.pack()
