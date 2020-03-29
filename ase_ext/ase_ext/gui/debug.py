import sys

import gtk


class Debug(gtk.Window):
    def __init__(self, gui):
        self.gui = gui
        gtk.Window.__init__(self)
        self.set_title('Debug')
        entry = gtk.Entry(200)
        self.add(entry)
        entry.connect('activate', self.enter, entry)
        entry.show()
        self.show()

    def enter(self, widget, entry):
        g = self.gui
        print(eval(entry.get_text()), file=sys.stderr)
    
