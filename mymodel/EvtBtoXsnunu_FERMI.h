/**************************************************************************
 * basf2 (Belle II Analysis Software Framework)                           *
 * Author: The Belle II Collaboration                                     *
 *                                                                        *
 * See git log for contributors and copyright holders.                    *
 * This file is licensed under LGPL-3.0, see LICENSE.md.                  *
 **************************************************************************/
#pragma once

#include "EvtGenBase/EvtDecayIncoherent.hh"
#include "EvtGenBase/EvtParticle.hh"

namespace Belle2 {
  class EvtBtoXsnunu_FERMI : public  EvtDecayIncoherent {

  public:

    EvtBtoXsnunu_FERMI() {}
    virtual ~EvtBtoXsnunu_FERMI();

    std::string getName();

    EvtDecayBase* clone();

    void initProbMax();

    void init();

    void decay(EvtParticle *p);

    double dGdsbProb(double _sb);
    double FermiMomentum(double pf);
    double FermiMomentumProb(double pb, double pf);

  private:
      double         _dGdsbProbMax;
      double         _mb;          // b-quark mass for fermi motion
      double         _ms;          // s-quark mass for fermi motion
      double         _mq;          // spectator quark mass for fermi motion
      double         _pf;          // Fermi motion parameter for fermi motion
      double         _mxmin;

  };

} // Belle 2 Namespace
