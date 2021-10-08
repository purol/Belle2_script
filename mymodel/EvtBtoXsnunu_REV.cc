// modified from EvtBtoXsll.cpp
// updated: 2021-10-05

#include <generators/evtgen/EvtGenModelRegister.h>

#include <stdlib.h>
#include <cmath>
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtId.hh"

#include "generators/evtgen/models/EvtBtoXsnunu_REV.h"

using std::endl;

namespace Belle2 {

  /** register the model in EvtGen */
  B2_EVTGEN_REGISTER_MODEL(EvtBtoXsnunu_REV);

  EvtBtoXsnunu_REV::~EvtBtoXsnunu_REV() {}

  std::string EvtBtoXsnunu_REV::getName()
  {
    return "BTOXSNUNU_REV";
  }

  EvtDecayBase* EvtBtoXsnunu_REV::clone()
  {
    return new EvtBtoXsnunu_REV;
  }

  void EvtBtoXsnunu_REV::decay(EvtParticle* p)
  {
      static EvtId NUE = EvtPDL::getId("nu_e");
      static EvtId NUM = EvtPDL::getId("nu_mu");
      static EvtId NUT = EvtPDL::getId("nu_tau");
      static EvtId NUEB = EvtPDL::getId("anti-nu_e");
      static EvtId NUMB = EvtPDL::getId("anti-nu_mu");
      static EvtId NUTB = EvtPDL::getId("anti-nu_tau");

      p->initializePhaseSpace(getNDaug(), getDaugs());

      double m_B = p->mass();

      EvtParticle* Xs, * neutrino1, * neutrino2;
      Xs = p->getDaug(0);
      neutrino1 = p->getDaug(1);
      neutrino2 = p->getDaug(2);
      EvtVector4R momnu1 = neutrino1->getP4();
      EvtVector4R momnu2 = neutrino2->getP4();
      EvtVector4R momXs = Xs->getP4();

      EvtVector4R q = momnu1 + momnu2;
      double q2 = q.mass2();

      double mb = 4.68;
      double ms = 0.1;
      double mstilda = ms / mb;

      double sb = q2 / (mb * mb);
      double mstilda2 = mstilda * mstilda;

      double lambda = 1 + mstilda2 * mstilda2 + sb * sb - 2 * (mstilda2 + sb + mstilda2 * sb);

      double prob = sqrt(lambda) * (3 * sb * (1 + mstilda2 - sb) + lambda);

      setProb(prob);

      return;

  }


  void EvtBtoXsnunu_REV::initProbMax()
  {
      double mb = 4.68;
      double ms = 0.1;
      double mstilda = ms / mb;
      double mstilda2 = mstilda * mstilda;

      int nsteps = 100;
      double sb = 0.0;
      double sbmin = 0;
      double sbmax = (1 - mstilda) * (1 - mstilda);
      double probMax = -10000.0;
      double sbProbMax = -10.0;

      for (int i = 0; i < nsteps; i++)
      {
          sb = sbmin + (i + 0.0005) * (sbmax - sbmin) / (double)nsteps;
          double lambda = 1 + mstilda2 * mstilda2 + sb * sb - 2 * (mstilda2 + sb + mstilda2 * sb);
          double prob = sqrt(lambda) * (3 * sb * (1 + mstilda2 - sb) + lambda);
          if (prob > probMax)
          {
              sbProbMax = sb;
              probMax = prob;
          }
      }

      if (verbose()) {
          std::cout << "dGdsbProbMax = " << probMax << " for sb = " << sbProbMax << std::endl;
      }

      setProbMax(probMax);
  }


  void EvtBtoXsnunu_REV::init()
  {

      // check that there are no arguments

      checkNArg(0);

      checkNDaug(3);

      // Check that the two leptons are the same type

      EvtId lepton1type = getDaug(1);
      EvtId lepton2type = getDaug(2);

      int etyp = 0;
      int mutyp = 0;
      int tautyp = 0;
      if (lepton1type == EvtPDL::getId("anti-nu_e") ||
          lepton1type == EvtPDL::getId("nu_e")) {
          etyp++;
      }
      if (lepton2type == EvtPDL::getId("anti-nu_e") ||
          lepton2type == EvtPDL::getId("nu_e")) {
          etyp++;
      }
      if (lepton1type == EvtPDL::getId("anti-nu_mu") ||
          lepton1type == EvtPDL::getId("nu_mu")) {
          mutyp++;
      }
      if (lepton2type == EvtPDL::getId("anti-nu_mu") ||
          lepton2type == EvtPDL::getId("nu_mu")) {
          mutyp++;
      }
      if (lepton1type == EvtPDL::getId("anti-nu_tau") ||
          lepton1type == EvtPDL::getId("nu_tau")) {
          tautyp++;
      }
      if (lepton2type == EvtPDL::getId("anti-nu_tau") ||
          lepton2type == EvtPDL::getId("nu_tau")) {
          tautyp++;
      }

      if (etyp != 2 && mutyp != 2 && tautyp != 2) {

          std::cout << "Expect two neutrinos of the same type in EvtBtoXsnunu.cc\n";
          ::abort();
      }

      // Check that the second and third entries are leptons with positive
      // and negative charge, respectively

      int lpos = 0;
      int lneg = 0;
      if (lepton1type == EvtPDL::getId("anti-nu_e") ||
          lepton1type == EvtPDL::getId("anti-nu_mu") ||
          lepton1type == EvtPDL::getId("anti-nu_tau")) {
          lpos++;
      }
      if (lepton2type == EvtPDL::getId("nu_e") ||
          lepton2type == EvtPDL::getId("nu_mu") ||
          lepton2type == EvtPDL::getId("nu_tau")) {
          lneg++;
      }

      if (lpos != 1 || lneg != 1) {

          std::cout << "Expect 2nd and 3rd particles to be anti-particle and particle of neutrinos in EvtBtoXsnunu.cc\n";
          ::abort();
      }
  }

}

