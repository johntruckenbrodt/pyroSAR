##############################################################
# graphical Interface for querying file names
# module of software gammaGUI
# John Truckenbrodt 2014-2015
##############################################################

from Tkinter import *
from tkFileDialog import askopenfilename, askdirectory, asksaveasfilename

from auxiliary import Environment


# FILE DIALOG WITH PATH DISPLAY AND BROWSE BUTTON
class FileQuery(Frame):
    def __init__(self, parent, name, option):
        Frame.__init__(self, parent, Environment.bcolor)
        self.pack()
        self.option = option
        self.lab = Label(self, Environment.label_ops, text=name, width=20)
        self.lab.pack(side=LEFT)
        self.file = StringVar()
        self.filename = StringVar()
        self.dirname = StringVar()
        self.ent = Entry(self, width=40, textvariable=self.file, highlightbackground="white")
        self.ent.pack(side=LEFT)
        self.btn = Button(self, Environment.button_ops, text="Browse", command=self.browse_button)
        self.btn.pack(side=LEFT)

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