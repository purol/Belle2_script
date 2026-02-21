/**************************************************************************
 * basf2 (Belle II Analysis Software Framework)                           *
 * Author: The Belle II Collaboration                                     *
 *                                                                        *
 * See git log for contributors and copyright holders.                    *
 * This file is licensed under LGPL-3.0, see LICENSE.md.                  *
 **************************************************************************/
#pragma once

#include "EvtGenBase/EvtDecayProb.hh"
#include "EvtGenBase/EvtParticle.hh"

namespace Belle2 {
  class EvtBtoXsnunu_REV: public  EvtDecayProb {

  public:

    EvtBtoXsnunu_REV() {}
    virtual ~EvtBtoXsnunu_REV();

    std::string getName();

    EvtDecayBase* clone();

    void initProbMax();

    void init();

    void decay(EvtParticle *p);

  };

} // Belle 2 Namespace
