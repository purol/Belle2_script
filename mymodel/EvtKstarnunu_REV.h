//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 1998      Caltech, UCSB
//
// Module: EvtGen/EvtKstarnunu.hh
//
// Description:
//
// Modification history:
//
//    DJL/RYD     August 11, 1998         Module created
//
//------------------------------------------------------------------------
// revised version of EvtKstarnunu
// updated 2021-10-05

#pragma once

#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtDecayAmp.hh"

class EvtParticle;

namespace Belle2 {
	class EvtKstarnunu_REV :public  EvtDecayAmp {

	public:

		EvtKstarnunu_REV() {}
		virtual ~EvtKstarnunu_REV();

		std::string getName();
		EvtDecayBase* clone();

		void init();

		void decay(EvtParticle* p);

	};

} // Belle 2 Namespace
