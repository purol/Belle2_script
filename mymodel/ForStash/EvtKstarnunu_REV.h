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
  class EvtKstarnunu_REV : public  EvtDecayAmp {

    /** The evtgen model to produce B-> Kstar nu nubar decay sample.
     * From factors are based on [arXiv:1503.05534v3].
     * 
     * You can use the model as follows:
     * 
     * Decay MyB+
     * 1.0 K*+ nu_e anti-nu_e                    PHOTOS KSTARNUNU_REV 0.38 -1.17 2.42 0.3 0.39 1.19 0.27 0.53 0.48;
     * Enddecay
     * 
     * First three arguments are z expansion coefficients for V0.
     * Next three arguments are z expansion coefficients for A1.
     * Next three arguments are z expansion coefficients for A12
     */

  public:

    /**
     * Constructor.
     */
    EvtKstarnunu_REV() {}

    /**
     * Destructor.
     */
    virtual ~EvtKstarnunu_REV();

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
      /** z expansion coefficient alpha0_v0    */
      double m_alpha0_v0{ 0.38 };

      /** z expansion coefficient alpha1_v0    */
      double m_alpha1_v0{ -1.17 };

      /** z expansion coefficient alpha2_v0    */
      double m_alpha2_v0{ 2.42 };

      /** z expansion coefficient alpha0_A1    */
      double m_alpha0_A1{ 0.3 };

      /** z expansion coefficient alpha1_A1    */
      double m_alpha1_A1{ 0.39 };

      /** z expansion coefficient alpha2_A1    */
      double m_alpha2_A1{ 1.19 };

      /** z expansion coefficient alpha0_A12    */
      double m_alpha0_A12{ 0.27 };

      /** z expansion coefficient alpha1_A12    */
      double m_alpha1_A12{ 0.53 };

      /** z expansion coefficient alpha2_A12    */
      double m_alpha2_A12{ 0.48 };

  };

} // Belle 2 Namespace
