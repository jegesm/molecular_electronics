# encoding: utf-8

"Module for calculating energies and forces."

import gtk
from ase_ext.gui.simulation import Simulation
from ase_ext.gui.widgets import oops, pack

class OutputFieldMixin:
    def makeoutputfield(self, box, label="Output:"):
        frame = gtk.Frame(label)
        if box is not None:
            box.pack_start(frame, True, True, 0)
        box2 = gtk.VBox()
        frame.add(box2)
        #pack(box, [gtk.Label("Output:")])
        scrwin = gtk.ScrolledWindow()
        scrwin.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)
        self.output = gtk.TextBuffer()
        txtview = gtk.TextView(self.output)
        txtview.set_editable(False)
        scrwin.add(txtview)
        scrwin.show_all()
        box2.pack_start(scrwin, True, True, 0)
        self.savebutton = gtk.Button(stock=gtk.STOCK_SAVE)
        self.savebutton.connect('clicked', self.saveoutput)
        self.savebutton.set_sensitive(False)
        pack(box2, [self.savebutton])
        box2.show()
        frame.show()
        return frame
    
    def activate_output(self):
        self.savebutton.set_sensitive(True)

    def saveoutput(self, widget):
        chooser = gtk.FileChooserDialog(
            'Save output', None, gtk.FILE_CHOOSER_ACTION_SAVE,
            (gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL,
             gtk.STOCK_SAVE, gtk.RESPONSE_OK))
        ok = chooser.run()
        if ok == gtk.RESPONSE_OK:
            filename = chooser.get_filename()
            txt = self.output.get_text(self.output.get_start_iter(),
                                       self.output.get_end_iter())
            f = open(filename, "w")
            f.write(txt)
            f.close()
        chooser.destroy()
        
class EnergyForces(Simulation, OutputFieldMixin):
    def __init__(self, gui):
        Simulation.__init__(self, gui)

        self.set_default_size(-1, 300)
        vbox = gtk.VBox()
        self.packtext(vbox,
                      "Calculate potential energy and the force on all atoms")
        self.packimageselection(vbox)
        pack(vbox, gtk.Label(""))
        self.forces = gtk.CheckButton("Write forces on the atoms")
        self.forces.set_active(True)
        pack(vbox, [self.forces])
        pack(vbox, [gtk.Label("")])
        self.makeoutputfield(vbox)
        pack(vbox, gtk.Label(""))
        self.makebutbox(vbox)
        vbox.show()
        self.add(vbox)
        self.show()
        self.gui.register_vulnerable(self)
        
    def run(self, *args):
        if not self.setup_atoms():
            return
        self.begin()
        e = self.atoms.get_potential_energy()
        txt = "Potential Energy:\n"
        txt += "  %8.3f eV\n\n" % (e,)
        if self.forces.get_active():
            txt +="Forces:\n"
            forces = self.atoms.get_forces()
            for f in forces:
                txt += "  %8.3f, %8.3f, %8.3f eV/Å\n" % tuple(f)
        self.output.set_text(txt)
        self.activate_output()
        self.end()
                
    def notify_atoms_changed(self):
        "When atoms have changed, check for the number of images."
        self.setupimageselection()
