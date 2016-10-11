##############################################################
# graphical Interface for querying file names
# module of software gammaGUI
# John Truckenbrodt 2014-2015
##############################################################

import Tkinter
from tkFileDialog import askopenfilename, askdirectory, asksaveasfilename

from auxiliary import Environment


# FILE DIALOG WITH PATH DISPLAY AND BROWSE BUTTON
class FileQuery(Tkinter.Frame):
    def __init__(self, parent, name, option):
        Tkinter.Frame.__init__(self, parent, Environment.bcolor)
        self.pack()
        self.option = option
        self.lab = Tkinter.Label(self, Environment.label_ops, text=name, width=20)
        self.lab.pack(side=Tkinter.LEFT)
        self.file = Tkinter.StringVar()
        self.filename = Tkinter.StringVar()
        self.dirname = Tkinter.StringVar()
        self.ent = Tkinter.Entry(self, width=40, textvariable=self.file, highlightbackground="white")
        self.ent.pack(side=Tkinter.LEFT)
        self.btn = Tkinter.Button(self, Environment.button_ops, text="Browse", command=self.browse_button)
        self.btn.pack(side=Tkinter.LEFT)

    # create button for browsing folders/files
    def browse_button(self):
        if self.option == 1:
            self.filename = askopenfilename(parent=self, initialdir=Environment.workdir.get(), title="browse file")
            self.ent.insert(0, self.filename)
        elif self.option == 2:
            self.dirname = askdirectory(parent=self, initialdir=Environment.workdir.get(), title="browse folder")
            if self.dirname:
                self.file.set(self.dirname)
        elif self.option == 3:
            self.filename = asksaveasfilename(parent=self, initialdir=Environment.workdir.get(), title="browse file")
            self.ent.insert(0, self.filename)
