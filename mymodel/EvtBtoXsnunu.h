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
//#include "EvtGenModels/EvtBtoXsllUtil.hh"

class EvtBtoXsllUtil;

namespace Belle2 {
  class EvtBtoXsnunu: public  EvtDecayIncoherent  {

  public:

    EvtBtoXsnunu(): _calcprob(0) {}
    virtual ~EvtBtoXsnunu();

    std::string getName();

    EvtDecayBase* clone();

    void initProbMax();

    void init();

    void decay(EvtParticle *p);

  private:
    EvtBtoXsllUtil* _calcprob;
    double _dGdsProbMax;
    double _dGdsdupProbMax;
    double  _mb;
    double _ms;
    double _mq;
    double _pf;
    double _mxmin;

  };

} // Belle 2 Namespace
