'''
Prompt and Fake Rate methods as defined in AN-2010/261.
This method computes the number of prompt and fake events for a given
number of leptons with associated efficiencies and number of events
passing loose and tight selections.
For now, only the 3-lepton case is coded. The generic case will be added
later.

Author: Devin N. Taylor, UW-Madison
'''

class PromptFakeRate(object):
    '''A class to calculated the fake rate.'''
    def __init__(self,numLeptons,**kwargs):
       '''Initialize the FakeRate object given a number of leptons'''
       # get arguments
       
       # setup
       self.nl = numLeptons

    def setPromptRate(self,*promptRates):
        '''
        Setup the prompt rates.
        This is the probability of a prompt lepton that passes the loose
        selection to also pass the tight selection.
        '''
        if len(promptRates)!=self.nl:
            print "Error. Number of efficiencies must match number of leptons"
        else:
            self.p = promptRates

    def setFakeRate(self,*fakeRates):
        '''
        Setup the fake rates.
        This is the probability of a fake lepton that passes the loose selection
        to also pass the tight selection.
        '''
        if len(fakeRates)!=self.nl:
            print "Error. Number of efficiencies must match number of leptons"
        else:
            self.f = fakeRates

    def setEventNumbers(self,eventMap):
        '''
        Define the event map of tight and loose leptons.
        The structure is:
            eventMap = {
                'TTT' : number of all pass tight,
                'TTL' : number of first two pass tight and third fails,
                'TLT' : number of first and third pass tight and second fails,
                ...
            }
        The number of characters in the keys must match the number of leptons
        in the setup. 'T' means 'pass loose and tight' and 'L means 'pass loose
        and fail tight.'
        '''
        goodEventMap = True
        for key in eventMap:
            if len(key) != self.nl:
                goodEventMap = False
            for c in key:
                if c not in ['T','L']:
                    goodEventMap = False
        if not goodEventMap:
            print "Error. Keys must be strings with the same length as the number of leptons and include only 'T' and 'L'."
        else:
            self.eventMap = eventMap

    def getSignalEvents(self):
        '''
        Returns the number of signal events given the already initialized prompt and fake rates and event numbers.
        '''
        if not hasattr(self,'p'):
            print 'Error: First initialize the prompt rates.'
            return -1
        elif not hasattr(self,'f'):
            print 'Error: First initialize the fake rates.'
            return -1
        elif not hasattr(self,'eventMap'):
            print 'Error: First initialize the event numbers'
            return -1
        else:
            return self.__calculateSignalEvents()

    def __calculateSignalEvents(self):
        '''
        The actual calculation. Only verified for 3 leptons at the moment.
        '''
        if self.nl != 3:
            print 'Error: only 3 leptons implemented.'
            return -1
        p = self.p
        f = self.f
        coeff = 1.
        for i in range(self.nl):
            coeff *= 1./(p[i]-f[i])
        val = 0
        for k,n in self.eventMap.iteritems():
            tempVal = 1.
            for i in range(self.nl):
                c = k[i]
                if c == 'T':
                    tempVal *= (1-f[i])
                if c == 'L':
                    tempVal *= (-f[i])
            val += tempVal
        return coeff * val


