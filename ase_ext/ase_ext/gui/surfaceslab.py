# encoding: utf-8
"""surfaceslab.py - Window for setting up surfaces
"""

import gtk
from ase_ext.gui.widgets import pack, cancel_apply_ok, oops
from ase_ext.gui.pybutton import PyButton
from ase_ext.gui.setupwindow import SetupWindow
import ase_ext.lattice.surface as _surf
import ase
import numpy as np

introtext = """\
  Use this dialog to create surface slabs.  Select the element by
writing the chemical symbol or the atomic number in the box.  Then
select the desired surface structure.  Note that some structures can
be created with an othogonal or a non-orthogonal unit cell, in these
cases the non-orthogonal unit cell will contain fewer atoms.

  If the structure matches the experimental crystal structure, you can
look up the lattice constant, otherwise you have to specify it
yourself."""

# Name, structure, orthogonal, support-nonorthogonal, function
surfaces = [('FCC(100)', 'fcc', True, False, _surf.fcc100),
            ('FCC(110)', 'fcc', True, False, _surf.fcc110),
            ('FCC(111) non-orthogonal', 'fcc', False, True, _surf.fcc111),
            ('FCC(111) orthogonal', 'fcc', True, True, _surf.fcc111),
            ('BCC(100)', 'bcc', True, False, _surf.bcc100),
            ('BCC(110) non-orthogonal', 'bcc', False, True, _surf.bcc110),
            ('BCC(110) orthogonal', 'bcc', True, True, _surf.bcc110),
            ('BCC(111) non-orthogonal', 'bcc', False, True, _surf.bcc111),
            ('BCC(111) orthogonal', 'bcc', True, True, _surf.bcc111),
            ('HCP(0001) non-orthogonal', 'hcp', False, True, _surf.hcp0001),
            ('HCP(0001) orthogonal', 'hcp', True, True, _surf.hcp0001),
            ]

py_template = """
from ase_ext.lattice.surface import %(func)s

atoms = %(func)s(symbol='%(symbol)s', size=%(size)s,
    a=%(a).3f, vacuum=%(vacuum).3f%(orthoarg)s)
"""

class SetupSurfaceSlab(SetupWindow):
    """Window for setting up a surface."""
    def __init__(self, gui):
        SetupWindow.__init__(self)
        self.set_title("Surface")
        self.atoms = None

        vbox = gtk.VBox()

        # Intoductory text
        self.packtext(vbox, introtext)
             
        # Choose the element
        label = gtk.Label("Element: ")
        element = gtk.Entry(max=3)
        self.element = element
        self.elementinfo = gtk.Label("")
        pack(vbox, [label, element, self.elementinfo])
        self.element.connect('activate', self.update)
        self.legal_element = False
        
        # Choose the surface structure
        label = gtk.Label("Structure: ")
        self.structchoice = gtk.combo_box_new_text()
        self.surfinfo = {}
        for s in surfaces:
            assert len(s) == 5
            self.structchoice.append_text(s[0])
            self.surfinfo[s[0]] = s
        pack(vbox, [label, self.structchoice])
        self.structchoice.connect('changed', self.update)

        # Choose the lattice constant
        tbl = gtk.Table(2, 3)
        label = gtk.Label("Lattice constant: ")
        tbl.attach(label, 0, 1, 0, 1)
        vbox2 = gtk.VBox()          # For the non-HCP stuff
        self.vbox_hcp = gtk.VBox()  # For the HCP stuff.
        self.lattice_const = gtk.Adjustment(3.0, 0.0, 1000.0, 0.01)
        lattice_box = gtk.SpinButton(self.lattice_const, 10.0, 3)
        lattice_box.numeric = True
        pack(vbox2, [gtk.Label("a:"), lattice_box, gtk.Label("Å")])
        tbl.attach(vbox2, 1, 2, 0, 1)
        lattice_button = gtk.Button("Get from database")
        tbl.attach(lattice_button, 2, 3, 0, 1)
        # HCP stuff
        self.hcp_ideal = (8.0/3)**(1.0/3)
        self.lattice_const_c = gtk.Adjustment(self.lattice_const.value * self.hcp_ideal,
                                              0.0, 1000.0, 0.01)
        lattice_box_c = gtk.SpinButton(self.lattice_const_c, 10.0, 3)
        lattice_box_c.numeric = True
        pack(self.vbox_hcp, [gtk.Label("c:"), lattice_box_c, gtk.Label("Å")])
        self.hcp_c_over_a_format = "c/a: %.3f (%.1f %% of ideal)"
        self.hcp_c_over_a_label = gtk.Label(self.hcp_c_over_a_format % (self.hcp_ideal,
                                                                        100.0))
        pack(self.vbox_hcp, [self.hcp_c_over_a_label])
        tbl.attach(self.vbox_hcp, 1, 2, 1, 2)
        tbl.show_all()
        pack(vbox, [tbl])
        self.lattice_const.connect('value-changed', self.update)
        self.lattice_const_c.connect('value-changed', self.update)
        lattice_button.connect('clicked', self.get_lattice_const)
        pack(vbox, gtk.Label(""))

        # System size
        self.size = [gtk.Adjustment(1, 1, 100, 1) for i in range(3)]
        buttons = [gtk.SpinButton(s, 0, 0) for s in self.size]
        self.vacuum = gtk.Adjustment(10.0, 0, 100.0, 0.1)
        vacuum_box = gtk.SpinButton(self.vacuum, 0.0, 1)
        pack(vbox, [gtk.Label("Size: \tx: "), buttons[0],
                    gtk.Label(" unit cells")])
        pack(vbox, [gtk.Label("\t\ty: "), buttons[1],
                    gtk.Label(" unit cells")])
        pack(vbox, [gtk.Label("      \t\tz: "), buttons[2],
                    gtk.Label(" layers,  "),
                    vacuum_box, gtk.Label(" Å vacuum")])
        self.nosize = "\t\tNo size information yet."
        self.sizelabel = gtk.Label(self.nosize)
        pack(vbox, [self.sizelabel])
        for s in self.size:
            s.connect('value-changed', self.update)
        self.vacuum.connect('value-changed', self.update)
        pack(vbox, gtk.Label(""))

        # Buttons
        self.pybut = PyButton("Creating a surface slab.")
        self.pybut.connect('clicked', self.update)
        buts = cancel_apply_ok(cancel=lambda widget: self.destroy(),
                               apply=self.apply,
                               ok=self.ok)
        pack(vbox, [self.pybut, buts], end=True, bottom=True)
        
        self.add(vbox)
        vbox.show()
        self.show()
        self.gui = gui

        # Hide the HCP stuff to begin with.
        self.vbox_hcp.hide_all()

    # update_element inherited from SetupWindow

    def update(self, *args):
        "Called when something has changed."
        struct = self.structchoice.get_active_text()
        if struct:
            structinfo = self.surfinfo[struct]
            if structinfo[1] == 'hcp':
                self.vbox_hcp.show_all()
                ca = self.lattice_const_c.value / self.lattice_const.value
                self.hcp_c_over_a_label.set_text(self.hcp_c_over_a_format %
                                                 (ca, 100 * ca / self.hcp_ideal))
            else:
                self.vbox_hcp.hide_all()
        # Abort if element or structure is invalid
        if not (self.update_element() and struct):
            self.sizelabel.set_text(self.nosize)
            self.atoms = None
            self.pybut.python = None
            return False
        # Make the atoms
        assert self.legal_element
        kw = {}
        kw2 = {}
        if structinfo[3]:  # Support othogonal keyword?
            kw['orthogonal'] = structinfo[2]
            kw2['orthoarg'] = ', orthogonal='+str(kw['orthogonal'])
        else:
            kw2['orthoarg'] = ''
        kw2['func'] = structinfo[4].__name__
        kw['symbol'] = self.legal_element
        kw['size'] = [int(s.value) for s in self.size]
        kw['a'] = self.lattice_const.value
        kw['vacuum'] = self.vacuum.value
        # Now create the atoms
        try:
            self.atoms = structinfo[4](**kw)
        except ValueError as e:
            # The values were illegal - for example some size
            # constants must be even for some structures.
            self.pybut.python = None
            self.atoms = None
            self.sizelabel.set_text(str(e).replace(".  ", ".\n"))
            return False
        kw2.update(kw)
        self.pybut.python = py_template % kw2
        # Find the heights of the unit cell
        h = np.zeros(3)
        uc = self.atoms.get_cell()
        for i in range(3):
            norm = np.cross(uc[i-1], uc[i-2])
            norm /= np.sqrt(np.dot(norm, norm))
            h[i] = np.abs(np.dot(norm, uc[i]))
        natoms = len(self.atoms)
        txt = ("\t\t%.2f Å x %.2f Å x %.2f Å,  %i atoms."
               % (h[0], h[1], h[2], natoms))
        self.sizelabel.set_text(txt)
        return True
    
    def get_lattice_const(self, *args):
        if not self.update_element():
            oops("Invalid element.")
            return
        z = ase_ext.atomic_numbers[self.legal_element]
        ref = ase_ext.data.reference_states[z]
        surface = self.structchoice.get_active_text()
        if not surface:
            oops("No structure specified!")
            return
        struct = self.surfinfo[surface][1]
        if ref is None or ref['symmetry'].lower() != struct:
            oops(struct.upper() + " lattice constant unknown for "
                      + self.legal_element + ".")
            return
        a = ref['a']
        self.lattice_const.set_value(a)
        if struct == 'hcp':
            c = ref['c/a'] * a
            self.lattice_const_c.set_value(c)

    def apply(self, *args):
        self.update()
        if self.atoms is not None:
            self.gui.new_atoms(self.atoms)
            return True
        else:
            oops("No valid atoms.",
                 "You have not (yet) specified a consistent set of parameters.")
            return False

    def ok(self, *args):
        if self.apply():
            self.destroy()
            

        
