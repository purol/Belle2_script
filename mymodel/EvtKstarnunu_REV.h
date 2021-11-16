// revised version of EvtKstarnunu
// updated 2021-10-05
// form factor: arXiv:1503.05534v3

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
