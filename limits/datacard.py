_EPSILON = 1.0e-10

class Datacard(object):

    def __init__(self, name):
        self.card_name = name
        self.bkg = []
        self.syst = []
        self.observed = 0

    def add_sig(self, name, rate):
        self.signal = (name, rate)

    def set_observed(self, obs):
        self.observed = obs

    def add_bkg(self, name, rate):
        if rate == 0.0:
            rate = _EPSILON
        self.bkg.append((name, rate))

    def add_syst(self, name, syst_type, **chans):
        self.syst.append((name, syst_type, chans))

    def dump(self):
        imax = 1
        jmax = len(self.bkg)
        kmax = len(self.syst)

        out = ""
        out += "#%s\n" % self.card_name
        out += "imax %2i number of channels\n" % imax
        out += "jmax %2i number of backgrounds\n" % jmax
        out += "kmax %2i number of nuisance parameters\n" % kmax

        out += "-" * 20 + "\n"

        out += "bin 1\n"
        out += "observation %i\n" % self.observed

        out += "-" * 20 + "\n"

        fmt = "{:<25}" + "{:^11}" * (1 + jmax) + "\n"
        fmt_f = "{:<25}" + "{:^11.3e}" * (1 + jmax) + "\n"

        row = ["1" for i in xrange(1+jmax)]
        out += fmt.format("bin", *row)

        row = [self.signal[0]] + [x[0] for x in self.bkg]
        out += fmt.format("process", *row)

        row = [i for i in xrange(1+jmax)]
        out += fmt.format("process", *row)

        row = [self.signal[1]] + [x[1] for x in self.bkg]
        out += fmt_f.format("rate", *row)

        out += "-" * 20 + "\n"

        fmt = "{:<20}" + "{:<5}" + "{:^11}" * (1 + jmax) + "\n"

        for syst in self.syst:
            names = [self.signal[0]] + [x[0] for x in self.bkg]
            row = [syst[2][name] if name in syst[2] else '-' for name in names]
            out += fmt.format(syst[0], syst[1], *row)

        return out

