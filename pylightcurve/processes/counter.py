
import sys
import time
import datetime


class Counter:

    def __init__(self, counter, total_iterations, show_every=1):

        self.counter = counter
        self.current_iteration = 0
        self.total_iterations = int(total_iterations)
        self.start_time = time.time()
        self.show = 0
        self.show_every = int(show_every)

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

            sys.stdout.write('\r\033[K')
            sys.stdout.write('{0}{1}: {2}   time left: {3}   total time: {4}'.format(
                self.counter, '.' * (15 - len(self.counter)), out_of, time_left, total_time))
            sys.stdout.flush()

            self.show = 0

        if self.current_iteration == self.total_iterations and self.total_iterations > 1:

            print('')
