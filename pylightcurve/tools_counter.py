from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from ._1databases import *


def initialise_window(window, window_name, windows_to_hide, windows_to_close, exit_python, deactivate_close=False):

    def exit_command():

        for i in windows_to_close:
            i.destroy()

        for i in windows_to_hide:
            i.withdraw()

        if exit_python:
            os._exit(-1)

    if deactivate_close:

        def exit_command():
            pass

    window.wm_title(window_name)
    window.protocol('WM_DELETE_WINDOW', exit_command)

    window.withdraw()


def setup_window(window, objects, main_font=None, button_font=None, entries_bd=3):

    import sys
    if sys.version_info[0] > 2:
        from tkinter import Label
    else:
        from Tkinter import Label

    if button_font is None:
        button_font = ['times', 15, 'bold']

    if main_font is None:
        main_font = ['times', 15]

    for row in range(len(objects)):
        if len(objects[row]) == 0:
            label_empty = Label(window, text='')
            label_empty.grid(row=row, column=100)
        else:
            for obj in objects[row]:

                if obj[0].winfo_class() == 'Button':
                    obj[0].configure(font=button_font)
                elif obj[0].winfo_class() == 'Entry':
                    obj[0].configure(bd=entries_bd, font=main_font)
                elif obj[0].winfo_class() in ['Label', 'Radiobutton']:
                    obj[0].configure(font=main_font)

                if len(obj) == 4:
                    obj[0].grid(row=row, column=obj[1], columnspan=obj[2], rowspan=obj[3])
                elif len(obj) == 3:
                    obj[0].grid(row=row, column=obj[1], columnspan=obj[2])
                else:
                    obj[0].grid(row=row, column=obj[1])


def finalise_window(window, position=5, topmost=False):

    window.update_idletasks()

    if position == 1:
        x = 0
        y = 0

    elif position == 2:
        x = (window.winfo_screenwidth() - window.winfo_reqwidth()) / 2
        y = 0

    elif position == 3:
        x = window.winfo_screenwidth() - window.winfo_reqwidth()
        y = 0

    elif position == 4:
        x = 0
        y = (window.winfo_screenheight() - window.winfo_reqheight()) / 2

    elif position == 5:
        x = (window.winfo_screenwidth() - window.winfo_reqwidth()) / 2
        y = (window.winfo_screenheight() - window.winfo_reqheight()) / 2

    elif position == 6:
        x = window.winfo_screenwidth() - window.winfo_reqwidth()
        y = (window.winfo_screenheight() - window.winfo_reqheight()) / 2

    elif position == 7:
        x = 0
        y = window.winfo_screenheight() - window.winfo_reqheight()

    elif position == 8:
        x = (window.winfo_screenwidth() - window.winfo_reqwidth()) / 2
        y = window.winfo_screenheight() - window.winfo_reqheight()

    elif position == 9:
        x = window.winfo_screenwidth() - window.winfo_reqwidth()
        y = window.winfo_screenheight() - window.winfo_reqheight()

    else:
        x = 0
        y = 0

    window.geometry('+%d+%d' % (x, y))

    window.update_idletasks()

    window.lift()
    window.wm_attributes("-topmost", 1)
    if not topmost:
        window.after_idle(window.attributes, '-topmost', 0)

    window.deiconify()


class Counter():

    def __init__(self, counter, counter_window, total_iterations, show_every=1):

        self.counter = counter

        if self.counter is not False:
            if isinstance(counter, str):
                self.counter = counter
            else:
                self.counter = 'Counter'

        self.counter_window = counter_window

        if self.counter_window is not False:
            if isinstance(counter_window, str):
                self.counter_window = counter_window
            else:
                self.counter_window = 'Counter'

        self.current_iteration = 0
        self.total_iterations = int(total_iterations)
        self.start_time = time.time()
        self.show = 0
        self.show_every = int(show_every)

        if self.total_iterations == 1:
            self.show_every = 10

        if self.counter_window:

            import sys
            if sys.version_info[0] > 2:
                from tkinter import Tk, Label
            else:
                from Tkinter import Tk, Label

            self.root = Tk()

            self.label1 = Label(self.root, text=self.counter_window)
            self.label2 = Label(self.root, text='COMPLETE:')
            self.label3 = Label(self.root, text=' ')
            self.label4 = Label(self.root, text='TIME LEFT:')
            self.label5 = Label(self.root, text=' ')
            self.label6 = Label(self.root, text='TOTAL TIME:')
            self.label7 = Label(self.root, text=' ')

            initialise_window(self.root, self.counter_window, [], [self.root], False, deactivate_close=True)

            setup_window(self.root, [
                [[self.label1, 1, 2]],
                [[self.label2, 1], [self.label3, 2]],
                [[self.label4, 1], [self.label5, 2]],
                [[self.label6, 1], [self.label7, 2]],
            ])

    def update(self, message=''):

        self.current_iteration += 1
        self.show += 1.0 / self.show_every

        out_of = '{0}{1} / {2} {3}'.format(' ' * (len(str(self.total_iterations)) - len(str(self.current_iteration))),
                                       str(self.current_iteration), str(self.total_iterations), message)

        delta_time = time.time() - self.start_time

        time_left = str(datetime.timedelta(
            seconds=int((self.total_iterations - self.current_iteration) * delta_time / self.current_iteration)))

        total_time = str(datetime.timedelta(seconds=int(delta_time)))

        if int(self.show):

            if self.counter_window:
                self.label3.configure(text='{0}{1}{2}'.format('     ', out_of, '     '))
                self.label5.configure(text='{0}{1}{2}'.format('     ', time_left, '     '))
                self.label7.configure(text='{0}{1}{2}'.format('     ', total_time, '     '))
                self.root.update()

            if self.counter:
                sys.stdout.write('\r\033[K')
                sys.stdout.write('{0}{1}: {2}   time left: {3}   total time: {4}'.format(
                    self.counter, '.' * (15 - len(self.counter)), out_of, time_left, total_time))
                sys.stdout.flush()

            self.show = 0

        if self.current_iteration == 1:

            if self.counter_window:
                self.label3.configure(text='{0}{1}{2}'.format('     ', out_of, '     '))
                self.label5.configure(text='{0}{1}{2}'.format('     ', time_left, '     '))
                self.label7.configure(text='{0}{1}{2}'.format('     ', total_time, '     '))
                self.root.update()

                finalise_window(self.root, topmost=True)

                self.root.update()

        elif self.current_iteration == self.total_iterations and self.total_iterations > 1:

            if self.counter_window:
                self.root.destroy()

            if self.counter:
                print('')
