import pypdt
import operator

class Investigator:
    def __init__(self):
        self.FirstBDecayString = []
        self.SecondBDecayString = []
        self.FirstBFinalParticles = {}
        self.SecondBFinalParticles = {}
    def PutString(self, string, Type):
        if Type == "firstB":
            self.FirstBDecayString.append(string)
        elif Type == "secondB":
            self.SecondBDecayString.append(string)
        else:
            print("Invalid Type!")
            exit()
    def Clean(self):
        self.FirstBDecayString = []
        self.SecondBDecayString = []
        self.FirstBFinalParticles = {}
        self.SecondBFinalParticles = {}
    def Investigate(self, Type):
        if Type == "firstB":
            DecayString = self.FirstBDecayString
            FinalParticleList = self.FirstBFinalParticles
        elif Type == "secondB":
            DecayString = self.SecondBDecayString
            FinalParticleList = self.SecondBFinalParticles
        else:
            print("Invalid type!")
            exit()
        FirstPID = DecayString[0].replace(" ", "")
        if not ((FirstPID == "521\n") or (FirstPID == "-521\n") or (FirstPID == "511\n") or (FirstPID == "-511\n")):
            print("PID of first particle in decay string is not B meson!")
            exit()
        for idx, val in enumerate(DecayString):
            if idx == len(DecayString) - 1: # last particle
                self.AddInFinalParticles( val.replace(" ", ""), FinalParticleList )
            else:
                NWhiteThis = val.count(" ")
                NWhiteNext = val.count(" ")
                if NWhiteThis >= NWhiteNext: # do not decay into other particles
                    self.AddInFinalParticles( val.replace(" ", ""), FinalParticleList )
    def InvestigateAll(self):
        self.Investigate("firstB")
        self.Investigate("secondB")
    def AddInFinalParticles(self, PID, dictionary):
        if PID in dictionary:
            dictionary[PID] = dictionary[PID] +1
        else:
            dictionary[PID] = 1
    def GetParticles(self, DecayString, layer): # 4 space = 1 layer
        DecayParticles = ""
        
        # ignore charge conjugate
        Bmeson = DecayString[0].replace(" ", "")[0:-1]
        IsItConjugate = False
        if Bmeson[0] == "-":
            IsItConjugate = True
        if IsItConjugate:
            Bmeson = Bmeson[1:]
        DecayParticles = Bmeson + "=>"

        for String in DecayString:
            if ("    " * layer in String) and ("    " *layer + " " not in String):
                PID = String[4*layer:-1]
                if IsItConjugate:
                    if PID[0] == "-":
                        PID = PID[1:]
                    else:
                        PID = "-" + PID
                if PID == "-111": # exception pi0
                    PID = "111"
                elif PID == "-113": # exception rho(770)0
                    PID = "113"
                DecayParticles = DecayParticles + PID + "|"
        return DecayParticles
    def printBKG(self, ParticlePair):
        self.InvestigateAll()

        DecayParticles = self.GetParticles(self.FirstBDecayString, 2)
        if DecayParticles in ParticlePair:
            ParticlePair[DecayParticles] = ParticlePair[DecayParticles] + 1
        else:
            ParticlePair[DecayParticles] = 1

        DecayParticles = self.GetParticles(self.SecondBDecayString, 2)
        if DecayParticles in ParticlePair:
            ParticlePair[DecayParticles] = ParticlePair[DecayParticles] + 1
        else:
            ParticlePair[DecayParticles] = 1

with open("MIX_large_FBDT") as openfileobject:
    ParticlePair = {}
    MeetTagB = False
    MeetSigB = False
    MeetMCDecay = False
    UnderFirstB = False
    UnderSecondB = False
    Inv = Investigator()
    for line in openfileobject:
        if line == "\n":
            continue

        if (MeetTagB == False) and ((line.replace(" ", "") == "511\n") or (line.replace(" ", "") == "-511\n") or (line.replace(" ", "") == "521\n") or (line.replace(" ", "") == "-521\n")):
            MeetTagB = True
        elif (MeetTagB == True) and (MeetSigB == False) and (line.replace(" ", "") == "511\n" or line.replace(" ", "") == "-511\n" or line.replace(" ", "") == "521\n" or line.replace(" ", "") == "-521\n"):
            MeetSigB = True
        elif (MeetTagB == True) and (MeetSigB == True) and (line == "Monte Carlo Decay: \n"):
            MeetMCDecay = True
        elif (MeetMCDecay == True) and (UnderFirstB == False) and (line.replace(" ", "") == "511\n" or line.replace(" ", "") == "-511\n" or line.replace(" ", "") == "521\n" or line.replace(" ", "") == "-521\n"):
            UnderFirstB = True
        elif (MeetMCDecay == True) and (UnderFirstB == True) and (UnderSecondB == False) and (line.replace(" ", "") == "511\n" or line.replace(" ", "") == "-511\n" or line.replace(" ", "") == "521\n" or line.replace(" ", "") == "-521\n"):
            UnderSecondB = True
        elif line == "===============================\n":
            Inv.printBKG(ParticlePair)
            Inv.Clean()
            MeetTagB = False
            MeetSigB = False
            MeetMCDecay = False
            UnderFirstB = False
            UnderSecondB = False
        elif (line.replace(" ", "") == "511\n" or line.replace(" ", "") == "-511\n" or line.replace(" ", "") == "521\n" or line.replace(" ", "") == "-521\n"):
            print("unexpected B encounter!")
            exit()

        if (MeetTagB == True) and (MeetSigB == False):
            pass
        elif (MeetTagB == True) and (MeetSigB == True) and (MeetMCDecay == False):
            pass
        elif (MeetMCDecay == True) and (UnderFirstB == True) and (UnderSecondB == False):
            Inv.PutString(line, "firstB")
        elif (MeetMCDecay == True) and (UnderSecondB == True):
            Inv.PutString(line, "secondB")

print( sorted(ParticlePair.items(),key=operator.itemgetter(1), reverse = True) )

TotalN = 0
for Decay in ParticlePair:
    TotalN = TotalN + ParticlePair[Decay]
print(TotalN)
#f = open("CHG_large_FBDT", "r")

#f.readline()

#f.close()
