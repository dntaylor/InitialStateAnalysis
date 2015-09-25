_EPSILON = 1.0e-10

class Datacard(object):

    def __init__(self, name):
        self.card_name = name
        self.bkg = []
        self.syst = []
        self.signal = []
        self.observed = 0

    def add_sig(self, name, rate):
        self.signal.append((name, rate))

    def set_observed(self, obs):
        self.observed = obs

    def add_bkg(self, name, rate):
        if rate == 0.0:
            rate = _EPSILON
        self.bkg.append((name, rate))

    def add_syst(self, name, syst_type, **chans):
        self.syst.append((name, syst_type, chans))

    def dump(self):
        smax = len(self.signal)
        bmax = len(self.bkg)
        imax = 1
        jmax = smax
        kmax = len(self.syst)

        out = ""
        out += "#%s\n" % self.card_name
        out += "imax %2i number of channels\n" % imax
        out += "jmax %2i number of processes minus 1\n" % jmax
        out += "kmax %2i number of nuisance parameters\n" % kmax

        out += "-" * 30 + "\n"

        out += "bin 1\n"
        out += "observation %i\n" % self.observed

        out += "-" * 30 + "\n"

        fmt = "{:<40}" + "{:^11}" * (smax + bmax) + "\n"
        fmt_f = "{:<40}" + "{:^11.3e}" * (smax + bmax) + "\n"

        row = ["1" for i in xrange(smax+bmax)]
        out += fmt.format("bin", *row)

        row = [x[0] for x in self.signal + self.bkg]
        out += fmt.format("process", *row)

        row = [i-(smax-1) for i in xrange(smax+bmax)]
        out += fmt.format("process", *row)

        row = [x[1] for x in self.signal + self.bkg]
        out += fmt_f.format("rate", *row)

        out += "-" * 30 + "\n"

        fmt = "{:<30}" + "{:<10}" + "{:^11}" * (smax + bmax) + "\n"
        fmt_f = "{:<30}" + "{:<10}" + "{:^11.8e}" * (smax + bmax) + "\n"

        for syst in self.syst:
            names = [x[0] for x in self.signal + self.bkg]
            row = [syst[2][name] if name in syst[2] else '-' for name in names]
            #if 'gmN' in syst[1]:
            #    out += fmt_f.format(syst[0], syst[1], *row)
            #else:
            out += fmt.format(syst[0], syst[1], *row)

        return out

