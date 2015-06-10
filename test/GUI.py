__author__ = 'Alessandro Orchini'

# ! /usr/bin/python
"""Create GUI to setup parameters and run"""

from Tkinter import *

class Application():

    def __init__(self):
        """Initialize the application
        """
        #Default Flame Parameters
        beta = 6.0 # flame aspect ratio
        Mark = 0.02 # Markstein length
        axi = 1 # axisymmetric flag

        #Default Acoustic Parameters
        xf = 0.25 # flame position
        dT = 300. # finale temperature

        #Default Velocity Parameters
        K = 1.2 # Convection speed
        radComp = 1 # radial Component

        #Default Time Marching Parameters
        tf = 1.0 # final time

        #Default Transfer Function Parameters
        fMin = 25.0 # min freq
        fMax = 250.0 # max freq

        #Copy values in dictionaries
        #If I add a value here I need to update also the corresponding SetSomethingValue function
        self.flamePar = {'beta': beta, 'Mark': Mark, 'axi': axi}
        self.acousticPar = {'xf': xf, 'dT': dT}
        self.velocityPar = {'K': K, 'radComp': radComp}
        self.timePar = {'tf': tf}
        self.transferFunctionPar = {'fMin': fMin, 'fMax': fMax}

        root = Tk()
        menubar = self.createMenu(root)
        root.config(menu=menubar)

        print self.flamePar
        root.mainloop()
        root.destroy()


    def createMenu(self, root):
        """create a Menu
        """
        # Initialize the main menu
        menubar = Menu(root)
        filemenu = Menu(menubar, tearoff=0)

        # Add commands to this menu

        # File Menu
        filemenu.add_command(label="New", command=self.donothing) # This should reset everything
        filemenu.add_command(label="Open", command=self.donothing) # Load a saved state
        filemenu.add_command(label="Save", command=self.donothing) # Save the current state
        filemenu.add_separator() # Separator line. Estetic add on
        filemenu.add_command(label="Exit", command=root.quit) # Exit the program
        menubar.add_cascade(label="File", menu=filemenu ) # 	Creates a new hierarchical menu by associating a given menu to a parent menu

        # Edit Menu
        editmenu = Menu(menubar, tearoff=0)
        editmenu.add_command(label="Flame Parameters", command=self.setFlameParameters) # Set Flame parameters
        editmenu.add_command(label="Velocity Parameters", command=self.setVelocityParameters) # Set Velocity parameters
        editmenu.add_command(label="Acoustic Parameters", command=self.setAcousticParameters) # Set Acoustic parameters
        editmenu.add_command(label="Transfer Function Evaluation", command=self.setTransferFunctionParameters) # Set Acoustic parameters
        editmenu.add_command(label="Time Domain Simulations", command=self.setTimeMarchingParameters) # Set Acoustic parameters
        menubar.add_cascade(label="Edit", menu=editmenu)

        # Run Menu
        runmenu = Menu(menubar, tearoff=0)
        runmenu.add_command(label="Run Eigenvalue Problem", command=self.donothing) # Run eigp
        runmenu.add_command(label="Calculate Transfer Function", command=self.donothing) # TF evaluation
        runmenu.add_command(label="Run Time Domain Simulation", command=self.donothing) # Time domain calculations
        menubar.add_cascade(label="Edit", menu=runmenu)

        # Help Menu
        helpmenu = Menu(menubar, tearoff=0)
        helpmenu.add_command(label="Help Index", command=self.donothing)
        helpmenu.add_command(label="About...", command=self.donothing)
        menubar.add_cascade(label="Help", menu=helpmenu)

        return menubar

    def update_values(self, dict, entry):
        for key in dict:
            dict[key] = float(entry[key].get())


    def setFlameParameters(self):
        top = Toplevel()
        top.title("Edit Flame Parameters")

        # Loop over entries
        ind = 0
        entry = {}
        for key in self.flamePar:
            Label(top, text=key).grid(row = ind)
            entry[key] = Entry(top)
            entry[key].grid(row=ind, column=1)
            entry[key].insert(0, self.flamePar[key])
            ind = ind+1

        Button(top, text='Quit', command=top.quit).grid(row=ind, column=0, sticky=W, pady=4)
        Button(top, text='Update Values', command=(lambda a=self.flamePar, b=entry: self.update_values(a, b))).grid(row=ind, column=1, sticky=W, pady=4)
        mainloop()

        # Checks
        if(self.flamePar['axi']!=0 and self.flamePar['axi']!=1):
            print 'Error: the axi flag must be 0 or 1. Restored to default value '
            self.flamePar['axi'] = 1

        if(self.flamePar['Mark']<0):
            print 'Warning: negative Markstein length with this model can give non-physical results '

        if(self.flamePar['beta']<=0):
            print 'Error: beta must be positive. Restored to default value \n'
            self.flamePar['beta'] = 6.

        self.flamePar['axi'] = int(self.flamePar['axi'])
        print 'The updated Flame Variables are:'
        print self.flamePar
        print '\n'
        top.destroy()

    def setVelocityParameters(self):
        top = Toplevel()
        top.title("Edit Flame Parameters")

        # Loop over entries
        ind = 0
        entry = {}
        for key in self.velocityPar:
            Label(top, text=key).grid(row = ind)
            entry[key] = Entry(top)
            entry[key].grid(row=ind, column=1)
            entry[key].insert(0, self.velocityPar[key])
            ind = ind+1

        Button(top, text='Quit', command=top.quit).grid(row=ind, column=0, sticky=W, pady=4)
        Button(top, text='Update Values', command=(lambda a=self.velocityPar, b=entry: self.update_values(a, b))).grid(row=ind, column=1, sticky=W, pady=4)
        mainloop()

        # Checks
        if(self.velocityPar['radComp']!=0 and self.velocityPar['radComp']!=1):
            print 'Error: the radComp flag must be 0 or 1. Restored to default value '
            self.velocityPar['radComp'] = 1

        if(self.velocityPar['K']<=0):
            print 'Error: K must be positive. Restored to default value \n'
            self.velocityPar['K'] = 1.2

        self.velocityPar['radComp'] = int(self.velocityPar['radComp'])
        print 'The updated Velocity Variables are:'
        print self.velocityPar
        print '\n'
        top.destroy()

    def setAcousticParameters(self):
        top = Toplevel()
        top.title("Edit Flame Parameters")

        # Loop over entries
        ind = 0
        entry = {}
        for key in self.acousticPar:
            Label(top, text=key).grid(row = ind)
            entry[key] = Entry(top)
            entry[key].grid(row=ind, column=1)
            entry[key].insert(0, self.acousticPar[key])
            ind = ind+1

        Button(top, text='Quit', command=top.quit).grid(row=ind, column=0, sticky=W, pady=4)
        Button(top, text='Update Values', command=(lambda a=self.acousticPar, b=entry: self.update_values(a, b))).grid(row=ind, column=1, sticky=W, pady=4)
        mainloop()

        # Checks
        if(self.acousticPar['xf']<0 or self.acousticPar['xf']>1):
            print 'Error: xf must lie in [0,1] (nondimensional flame position). Restored to default value '
            self.acousticPar['xf'] = 0.25

        if(self.acousticPar['dT']<0):
            print 'Error: dT must be non-negative. Restored to default value \n'
            self.acousticPar['dT'] = 300.

        print 'The updated Acoustic Variables are:'
        print self.acousticPar
        print '\n'
        top.destroy()

    def setTimeMarchingParameters(self):
        top = Toplevel()
        top.title("Edit Flame Parameters")

        # Loop over entries
        ind = 0
        entry = {}
        for key in self.timePar:
            Label(top, text=key).grid(row = ind)
            entry[key] = Entry(top)
            entry[key].grid(row=ind, column=1)
            entry[key].insert(0, self.timePar[key])
            ind = ind+1

        Button(top, text='Quit', command=top.quit).grid(row=ind, column=0, sticky=W, pady=4)
        Button(top, text='Update Values', command=(lambda a=self.timePar, b=entry: self.update_values(a, b))).grid(row=ind, column=1, sticky=W, pady=4)
        mainloop()

        # Checks
        if(self.timePar['tf']<=0):
            print 'Error: tf must be positive. Restored to default value \n'
            self.timePar['tf'] = 1.

        print 'The updated Time Marching Variables are:'
        print self.timePar
        print '\n'
        top.destroy()

    def setTransferFunctionParameters(self):
        top = Toplevel()
        top.title("Edit Flame Parameters")

        # Loop over entries
        ind = 0
        entry = {}
        for key in self.transferFunctionPar:
            Label(top, text=key).grid(row = ind)
            entry[key] = Entry(top)
            entry[key].grid(row=ind, column=1)
            entry[key].insert(0, self.transferFunctionPar[key])
            ind = ind+1

        Button(top, text='Quit', command=top.quit).grid(row=ind, column=0, sticky=W, pady=4)
        Button(top, text='Update Values', command=(lambda a=self.transferFunctionPar, b=entry: self.update_values(a, b))).grid(row=ind, column=1, sticky=W, pady=4)
        mainloop()

        # Checks
        if(self.transferFunctionPar['fMin']<=0):
            print 'Error: fMin must be positive. Restored to default value \n'
            self.transferFunctionPar['fMin'] = 25.

        if(self.transferFunctionPar['fMax']<=0):
            print 'Error: fMax must be positive. Restored to default value \n'
            self.transferFunctionPar['fMax'] = 250.

        if(self.transferFunctionPar['fMin']>self.transferFunctionPar['fMax']):
            print 'Error: fMin < fMax. I will swap them'
            temp = self.transferFunctionPar['fMin']
            self.transferFunctionPar['fMin'] = self.transferFunctionPar['fMax']
            self.transferFunctionPar['fMax'] = temp

        if(self.transferFunctionPar['fMin'] == self.transferFunctionPar['fMax']):
            'Error: fMin = fMAx. This will note calculate a TF. For a single frequency analysis' \
            'use the eigenvalue problem solver. fMin and fMax are restored to their default values'
            self.transferFunctionPar['fMin'] = 25.
            self.transferFunctionPar['fMax'] = 250.

        print 'The updated Transfer Function Variables are:'
        print self.transferFunctionPar
        print '\n'
        top.destroy()

    def donothing(self):
       filewin = Toplevel()
       button = Button(filewin, text="Not yet implemented, sorry")
       button.pack()


Application()