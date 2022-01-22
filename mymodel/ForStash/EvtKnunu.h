/**************************************************************************
 * basf2 (Belle II Analysis Software Framework)                           *
 * Author: The Belle II Collaboration                                     *
 *                                                                        *
 * See git log for contributors and copyright holders.                    *
 * This file is licensed under LGPL-3.0, see LICENSE.md.                  *
 **************************************************************************/

#pragma once

#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtDecayAmp.hh"
#include <string>

class EvtParticle;

namespace Belle2 {
  class EvtKnunu : public  EvtDecayAmp {

    /** The evtgen model to produce B-> K nu nubar decay sample.
     * From factors are based on [arXiv:1409.4557v2].
     * 
     * You can use the model as follows:
     * 
     * Decay MyB+
     * 1.0 K+ nu_e anti-nu_e           PHOTOS KNUNU 0.432 -0.664 -1.2;
     * Enddecay
     * 
     * first, second, and third arguments are alpha0, alpha1, and alpha2 respectively for z expansion
     */

  public:

    /**
     * Constructor.
     */
    EvtKnunu() {}

    /**
     * Destructor.
     */
    virtual ~EvtKnunu();

    /**
     * The function which returns the name of the model.
     */
    std::string getName();

    /**
     * The function which makes a copy of the model.
     */
    EvtDecayBase* clone();

    /**
     * The function for an initialization.
     */
    void init();

    /**
     * The function to calculate a quark decay amplitude.
     */
    void decay(EvtParticle* p);

    /**
     * The function to sets a maximum probability.
     */
    void initProbMax();

  private:
      /** z expansion coefficient alpha0    */
      double m_alpha0{ 0.432 };

      /** z expansion coefficient alpha1    */
      double m_alpha1{ -0.664 };

      /** z expansion coefficient alpha2    */
      double m_alpha2{ -1.2 };

  };

} // Belle 2 Namespace
