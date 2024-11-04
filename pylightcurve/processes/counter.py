
import sys
import time
import datetime


class Counter:

    def __init__(self, counter, total_iterations, ignore=0, show_every=1, increment=1, notebook=False):

        self.counter = counter
        self.current_iteration = 0
        self.total_iterations = int(total_iterations)
        self.start_time = time.time()
        self.show = 0
        self.show_every = int(show_every)
        self.ignore = ignore
        self.max_message_len = 0
        self.increment = increment
        self.text = ''
        self.notebook = notebook

    def update(self, message=''):

        self.current_iteration += self.increment
        self.show += self.increment

        out_of = '{0}{1} / {2} {3}{4}'.format(
            ' ' * (len(str(self.total_iterations)) - len(str(self.current_iteration))),
            str(self.current_iteration), str(self.total_iterations), message, ' ' * (self.max_message_len - len(message)))

        elapsed_time_second = time.time() - self.start_time
        time_left_seconds = (self.total_iterations - self.current_iteration) * elapsed_time_second / self.current_iteration
        total_time_seconds = elapsed_time_second + time_left_seconds

        time_left = str(datetime.timedelta(seconds=int(time_left_seconds)))
        elapsed_time = str(datetime.timedelta(seconds=int(elapsed_time_second)))
        total_time = str(datetime.timedelta(seconds=int(total_time_seconds)))

        self.text = '{0}{1}: {2}, time left: {3}, time elapsed: {4}, total time: {5}'.format(
            self.counter, '.' * (15 - len(self.counter)), out_of, time_left, elapsed_time, total_time)

        if message != '':
            if len(message) > self.max_message_len:
                self.max_message_len = len(message)

        if self.ignore:
            if self.current_iteration == self.ignore:
                self.reset()
                self.ignore = 0

        if self.show == self.show_every and self.current_iteration > 0:

            if self.notebook:
                from IPython.display import clear_output, display
                clear_output(wait=True)
                print(self.text)
            else:
                sys.stdout.write('\r\033[K')
                sys.stdout.write(self.text)
                sys.stdout.flush()

            self.show = 0

        if self.current_iteration == self.total_iterations and self.total_iterations > 1:

            print('')

    def reset(self):

        self.current_iteration = 0
        self.start_time = time.time()
        self.show = 0
