// decay for B->K nu nu
// updated: 2021-10-05

#include <generators/evtgen/EvtGenModelRegister.h>

#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include <iostream>
#include <string>
#include <cmath>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtDiracSpinor.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtVector4C.hh"


#include "generators/evtgen/models/EvtKnunu.h"

using std::endl;

namespace Belle2 {

    /** register the model in EvtGen */
    B2_EVTGEN_REGISTER_MODEL(EvtKnunu);

    EvtKnunu::~EvtKnunu() { }

    std::string EvtKnunu::getName()
    {
        return "KNUNU";
    }

    EvtDecayBase* EvtKnunu::clone()
    {
        return new EvtKnunu;
    }

    void EvtKnunu::decay(EvtParticle* p)
    {
        static EvtId NUE = EvtPDL::getId("nu_e");
        static EvtId NUM = EvtPDL::getId("nu_mu");
        static EvtId NUT = EvtPDL::getId("nu_tau");
        static EvtId NUEB = EvtPDL::getId("anti-nu_e");
        static EvtId NUMB = EvtPDL::getId("anti-nu_mu");
        static EvtId NUTB = EvtPDL::getId("anti-nu_tau");

        p->initializePhaseSpace(getNDaug(), getDaugs());

        double m_b = p->mass();

        EvtParticle* meson, * neutrino1, * neutrino2;
        meson = p->getDaug(0);
        neutrino1 = p->getDaug(1);
        neutrino2 = p->getDaug(2);
        EvtVector4R momnu1 = neutrino1->getP4();
        EvtVector4R momnu2 = neutrino2->getP4();
        EvtVector4R momk = meson->getP4();

        EvtVector4R q = momnu1 + momnu2;
        double q2 = q.mass2();

        EvtVector4R p4b; p4b.set(m_b, 0., 0., 0.);  // Do calcs in mother rest frame

        double m_k = meson->mass();
        double mkhat = m_k / (m_b);
        double shat = q2 / (m_b * m_b);

        EvtVector4R phat = (p4b + momk) / m_b;
        EvtVector4R qhat = q / m_b;

        // calculate form factor
        // Buras, Andrzej J., et al. Journal of High Energy Physics 2015.2 (2015): 1-39.
        double alpha0 = 0.432;
        double alpha1 = -0.664;
        double alpha2 = -1.2;
        double mp = m_b + 0.046;
        double tp = (m_b + m_k) * (m_b + m_k);
        double tm = (m_b - m_k) * (m_b - m_k);
        double t0 = tp * (1 - sqrt(1 - tm / tp));
        double z = (sqrt(tp - q2) - sqrt(tp - t0)) / (sqrt(tp - q2) + sqrt(tp - t0));
        double fp = (1 / (1 - q2 / (mp * mp))) * (alpha0 + alpha1 * z + alpha2 * z * z + (-alpha1 + 2 * alpha2) * z * z * z / 3);
        double fm = -fp * (1 - mkhat * mkhat) / shat;

        EvtVector4C T1;

        EvtComplex aprime;
        aprime = fp;
        EvtComplex bprime;
        bprime = fm;

        T1 = aprime * phat + bprime * qhat;

        EvtVector4C l;

        if (getDaug(1) == NUE || getDaug(1) == NUM || getDaug(1) == NUT) {
            l = EvtLeptonVACurrent(neutrino1->spParentNeutrino(),
                neutrino2->spParentNeutrino());
        }
        if (getDaug(1) == NUEB || getDaug(1) == NUMB || getDaug(1) == NUTB) {
            l = EvtLeptonVACurrent(neutrino2->spParentNeutrino(),
                neutrino1->spParentNeutrino());
        }

        vertex(l * T1);

    }

    void EvtKnunu::init()
    {
        // check that there are 0 arguments
        checkNArg(0);
        checkNDaug(3);

        //We expect the parent to be a scalar 
        //and the daughters to be K neutrino netrino

        checkSpinParent(EvtSpinType::SCALAR);

        checkSpinDaughter(0, EvtSpinType::SCALAR);
        checkSpinDaughter(1, EvtSpinType::NEUTRINO);
        checkSpinDaughter(2, EvtSpinType::NEUTRINO);
    }

}

