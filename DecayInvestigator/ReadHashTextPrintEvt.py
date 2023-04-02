import pypdt
import operator
from collections import defaultdict

def default_value_for_new_key():
    return -1
 
MagicTable = defaultdict(default_value_for_new_key)

MagicTable['521=>-423|-13|14|'] = 0
MagicTable['521=>-421|-13|14|'] = 1
MagicTable['521=>-421|-11|12|'] = 2
MagicTable['521=>213|-421|'] = 3
MagicTable['521=>-423|-11|12|'] = 4
MagicTable['521=>-421|211|'] = 5
MagicTable['521=>-423|20213|'] = 6
MagicTable['521=>-413|111|211|211|'] = 7
MagicTable['521=>-423|211|'] = 8
MagicTable['521=>-421|-211|211|211|'] = 9
MagicTable['521=>-421|113|211|'] = 10
MagicTable['521=>-423|213|'] = 11
MagicTable['521=>-421|-15|16|'] = 12
MagicTable['521=>-423|-15|16|'] = 13
MagicTable['521=>-423|-211|211|211|111|'] = 14
MagicTable['521=>20213|-421|'] = 15
MagicTable['521=>-421|431|'] = 16
MagicTable['521=>321|311|-311|'] = 17
MagicTable['521=>-15|16|'] = 18
MagicTable['521=>-421|321|'] = 19
MagicTable['521=>-423|211|-211|211|'] = 20
MagicTable['521=>-421|321|-311|'] = 21
MagicTable['521=>433|-421|'] = 22
MagicTable['521=>-423|321|'] = 23
MagicTable['521=>-413|211|211|'] = 24
MagicTable['521=>443|321|'] = 25
MagicTable['521=>-423|431|'] = 26
MagicTable['511=>-413|-13|14|'] = 100
MagicTable['511=>-413|-11|12|'] = 101
MagicTable['511=>-411|-13|14|'] = 102
MagicTable['511=>-413|211|111|'] = 103
MagicTable['511=>-411|-11|12|'] = 104
MagicTable['511=>-413|20213|'] = 105
MagicTable['511=>213|-411|'] = 106
MagicTable['511=>213|-413|'] = 107
MagicTable['511=>-413|211|211|-211|111|'] = 108
MagicTable['511=>-413|113|211|'] = 109
MagicTable['511=>-413|211|'] = 110
MagicTable['511=>-411|-15|16|'] = 111
MagicTable['511=>-411|211|'] = 112
MagicTable['511=>-413|-15|16|'] = 113
MagicTable['511=>-411|-211|211|211|'] = 114
MagicTable['511=>20213|-411|'] = 115
MagicTable['511=>-413|431|'] = 116
MagicTable['511=>311|-311|311|'] = 117
MagicTable['511=>433|-413|'] = 118
MagicTable['511=>433|-411|'] = 119
MagicTable['511=>-421|-211|211|'] = 120
MagicTable['511=>-411|431|'] = 121
MagicTable['511=>-413|433|'] = 122

class Investigator:
    def __init__(self):
        self.FirstBDecayString = []
        self.SecondBDecayString = []
        self.FirstBFinalParticles = {}
        self.SecondBFinalParticles = {}
        self.AllString = []
    def PutString(self, string, Type):
        if Type == "firstB":
            self.FirstBDecayString.append(string)
        elif Type == "secondB":
            self.SecondBDecayString.append(string)
        elif Type == "all":
            self.AllString.append(string)
        else:
            print("Invalid Type!")
            exit()
    def Clean(self):
        self.FirstBDecayString = []
        self.SecondBDecayString = []
        self.FirstBFinalParticles = {}
        self.SecondBFinalParticles = {}
        self.AllString = []
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
                elif PID == "-443": # exception J/psi(1S)
                    PID = "443"
                DecayParticles = DecayParticles + PID + "|"
        return DecayParticles
    def printBKG(self, ParticlePair):
        self.InvestigateAll()

        DecayParticles = self.GetParticles(self.FirstBDecayString, 2)
        if DecayParticles in ParticlePair:
            ParticlePair[DecayParticles] = ParticlePair[DecayParticles] + 1
        else:
            ParticlePair[DecayParticles] = 1
        DMID_1 = MagicTable[DecayParticles]

        DecayParticles = self.GetParticles(self.SecondBDecayString, 2)
        if DecayParticles in ParticlePair:
            ParticlePair[DecayParticles] = ParticlePair[DecayParticles] + 1
        else:
            ParticlePair[DecayParticles] = 1
        DMID_2 = MagicTable[DecayParticles]

        print(self.AllString[0], end='') 
        print(self.AllString[1], end='')
        print(self.AllString[2], end='')
        print(self.AllString[3], end='')
        print(self.AllString[4], end='')
        print(DMID_1)
        print(DMID_2)

with open("MIX_test") as openfileobject:
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
        Inv.PutString(line, "all")

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
