import pypdt
import operator

class Investigator:
    def __init__(self):
        self.TagDecayString = []
        self.SigDecayString = []
        self.FirstBDecayString = []
        self.SecondBDecayString = []
        self.TagFinalParticles = {}
        self.FirstBFinalParticles = {}
        self.SecondBFinalParticles = {}
    def PutString(self, string, Type):
        if Type == "tag":
            self.TagDecayString.append(string)
        elif Type == "sig":
            self.SigDecayString.append(string)
        elif Type == "firstB":
            self.FirstBDecayString.append(string)
        elif Type == "secondB":
            self.SecondBDecayString.append(string)
        else:
            print("Invalid Type!")
            exit()
    def Clean(self):
        self.TagDecayString = []
        self.SigDecayString = []
        self.FirstBDecayString = []
        self.SecondBDecayString = []
        self.TagFinalParticles = {}
        self.FirstBFinalParticles = {}
        self.SecondBFinalParticles = {}
    def Investigate(self, Type):
        if Type == "tag":
            DecayString = self.TagDecayString
            FinalParticleList = self.TagFinalParticles
        elif Type == "firstB":
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
        self.Investigate("tag")
        self.Investigate("firstB")
        self.Investigate("secondB")
    def AddInFinalParticles(self, PID, dictionary):
        if PID in dictionary:
            dictionary[PID] = dictionary[PID] +1
        else:
            dictionary[PID] = 1
    def GetNParticles(self, dic):
        TotalN = 0
        for key, value in dic.items():
            TotalN = TotalN + value
        return TotalN

    def Evaluate(self):
        TagPIDSet = set(self.TagFinalParticles.keys())
        FirstBPIDSet = set(self.FirstBFinalParticles.keys())
        SecondBPIDSet = set(self.SecondBFinalParticles.keys())

        FirstBScore = 0
        SecondBScore = 0

        for PID in (TagPIDSet | FirstBPIDSet):
            TagN = 0
            FirstBN = 0
            if PID in self.TagFinalParticles:
                TagN = self.TagFinalParticles[PID]
            else:
                TagN = 0
            if PID in self.FirstBFinalParticles:
                FirstBN = self.FirstBFinalParticles[PID]
            else:
                FirstBN = 0
            FirstBScore = FirstBScore + abs(TagN - FirstBN)

        for PID in (TagPIDSet | SecondBPIDSet):
            TagN = 0
            SecondBN = 0
            if PID in self.TagFinalParticles:
                TagN = self.TagFinalParticles[PID]
            else:
                TagN = 0
            if PID in self.SecondBFinalParticles:
                SecondBN = self.SecondBFinalParticles[PID]
            else:
                SecondBN = 0
            SecondBScore = SecondBScore + abs(TagN - SecondBN)

        NFinalParticleFirstB = self.GetNParticles(self.FirstBFinalParticles)
        NFinalParticleSecondB = self.GetNParticles(self.SecondBFinalParticles)

        if FirstBScore <= SecondBScore:
            if (FirstBScore/NFinalParticleFirstB) < 0.5:
                return 1
        elif FirstBScore > SecondBScore:
            if (SecondBScore/NFinalParticleSecondB) < 0.5:
                return -1
        return 0
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
        result = self.Evaluate()
        if result == 1:
            DecayParticles = self.GetParticles(self.FirstBDecayString, 2)
        elif result == -1:
            DecayParticles = self.GetParticles(self.SecondBDecayString, 2)
        elif result == 0:
            DecayParticles = "others"
        else:
            print("unexpected result!")
            exit()
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
            Inv.PutString(line, "tag")
        elif (MeetTagB == True) and (MeetSigB == True) and (MeetMCDecay == False):
            Inv.PutString(line, "sig")
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
