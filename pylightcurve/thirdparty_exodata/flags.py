""" Contains the flag class, this is desgined to hold types of information about the target. These can be things
entered by the catalogue when assuming values such as 'Temperature calculated' or personal tags like 'Priority Target'

These are designed to be attached to a planet, system, star or binary class in .flags
"""

allowedFlags = ['Calculated Temperature', 'Estimated Mass', 'Calculated SMA', 'Fake', 'Estimated Distance', 'Calculated Period']
"UBVJIHKLMN"
allowedFlags += ['Estimated magU', 'Estimated magB', 'Estimated magV', 'Estimated magJ', 'Estimated magI',
                 'Estimated magH', 'Estimated magK', 'Estimated magL', 'Estimated magM', 'Estimated magN']


class Flags(object):  # or tags? or lists?

    def __init__(self):

        self.flags = set()

    def addFlag(self, flag):

        if flag in allowedFlags:
            self.flags.add(flag)
        else:
            raise InvalidFlag

    def removeFlag(self, flag):

        self.flags.remove(flag)

    def __repr__(self):

        return 'Flags({0})'.format(str(self.flags)[4:-1])

    def __iter__(self):

        return iter(self.flags)


class InvalidFlag(BaseException):
    pass
