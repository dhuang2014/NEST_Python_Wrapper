# LibNEST

This repository serves as the private LUX copy of the NEST software package.

This file should be kept up to date, detailing the differences between the public 
and LUX repositories. If changes are properly recorded, LibNEST development
can provide feedback on NEST improvement for future public releases. 

For full instructions on how to use the NEST software
(nearly identical usage), please see the file README.md.

## Tracking Internal LUX Development of NEST

In order to keep track of changes made to the private LUX repo, differences 
(see "diff" command) between latest LibNEST and NEST releases will be shown in this section.

Please add future diffs just below (to the top of this section).

### Diff: NEST v2.0.0 release and latest LibNEST (July 19, 2018)

Summary:
* Implementation of the correct detector classes "LUX_Run3" and "LUX_Run4" in "testNEST.cpp"
	(rather than the default "DetectorExample_XENON10" detector).
* Inclusion of the "LUXSpikeCountLookup.h" needed for LUX-based spike-counting in 
	"NEST.cpp", as well as the function "GetSpikeRun4" function which implements
	this spike-counting method.
* Also in "NEST.cpp", "GetS1" is modified to call the LUX-based spike-counting function. 
* nest repo has ".travis.yml" that is not needed by LUX.
* LibNEST contains a "Deprecated_Files" folder contains files from the previous
	LibNEST repo (for future reference).

```
Only in nest: .travis.yml
Only in LibNEST: Deprecated_Files

diff LibNEST/NEST.cpp nest/NEST.cpp
3d2
< #include "LUXSpikeCountLookup.h"
652,655c651,652
<     //newSpike = GetSpike(Nph, smearPos[0], smearPos[1], smearPos[2],
<     //                    driftVelocity, dV_mid, scintillation);
< 		newSpike = GetSpikeRun4 (nHits, smearPos[0], smearPos[1], smearPos[2],
< 												driftVelocity, dV_mid, type_num);
---
>     newSpike = GetSpike(Nph, smearPos[0], smearPos[1], smearPos[2],
>                         driftVelocity, dV_mid, scintillation);
1068,1201d1064
< vector<double> NESTcalc::GetSpikeRun4 ( int nHits, double dx, double dy, double dz,
<   double driftSpeed, double dS_mid, INTERACTION_TYPE species ) {
<
< 	// Output spike variables
< 	vector<double> newSpike(2);
<
< 	int numSpikesTop, numSpikesBot;
< 	int S1spike_raw;
< 	double S1spike_OLc;
< 	double numPHETop, numPHEBot;
< 	double S1Top, S1Bot;
< 	double spikeCountCorrectionTop, spikeCountCorrectionBot;
< 	double stRatio, tau;
< 	const double vuv1_1to1_2 = 0.915; //correction between different gain versions in LUX
< 	const double coinWind = 100.; //S1 concidence window in ns
< 	const double spEff_Run4 = 0.3;
<
< 	numPHETop = (double) BinomFluct(nHits, topArrayProbabilityCoeffs[0]*pow(dS_mid,2.) + topArrayProbabilityCoeffs[1]*dS_mid + topArrayProbabilityCoeffs[2]);
< 	numPHEBot = (double) nHits - numPHETop;
<
< 	// Do various checks on the number of PHE in the top and bottom arrays,
< 	// and obtain the RAW/UNCORRECTED number of spikes in top and bottom array
< 	if( int(numPHETop)-2 > 99 || int(numPHEBot)-2 > 99) {
< 		//cerr << "Number of photons in one of the arrays is too large for the spike count to work. " <<
< 		//      "Spike count variables will all be set to zero." << endl;
< 		numSpikesTop = 0;
< 		numSpikesBot = 0;
< 	}
< 	else {
< 		if (numPHETop < 2) numSpikesTop = numPHETop;
< 		else {
< 			numSpikesTop = (int) floor( RandomGen::rndm()->rand_gauss(spikes_S1[int(numPHETop)-2], sdev_spikes_S1[int(numPHETop)-2]) + 0.5 );
< 			if (numSpikesTop > numPHETop) numSpikesTop = numPHETop;
< 			if (numSpikesTop > 80) {
< 				//cerr << "Number of spikes in top array exceeds the allowed number. " <<
< 				//        "Spike counts will all be set to zero." << endl;
< 				numSpikesTop = 0;
< 				numSpikesBot = 0;
< 			}
< 		}
< 		if (numPHEBot < 2) numSpikesBot = numPHEBot;
< 		else {
< 			numSpikesBot = (int) floor( RandomGen::rndm()->rand_gauss(spikes_S1[int(numPHEBot)-2], sdev_spikes_S1[int(numPHEBot)-2]) + 0.5 );
< 			if (numSpikesBot > numPHEBot) numSpikesBot = numPHEBot;
< 			if (numSpikesBot > 80) {
< 				//cerr << "Number of spikes in bottom array exceeds the allowed number. " <<
<         //     	"Spike counts will all be set to zero." << endl;
< 				numSpikesTop = 0;
< 				numSpikesBot = 0;
< 			}
< 		}
< 	}
< 	// Get the total number of RAW/UNCORRECTED spikes
< 	S1spike_raw = numSpikesTop + numSpikesBot;
<
<
< 	// If there is a non-zero coincidence requirement, proceed
< 	if ( fdetector->get_coinLevel() ) {
< 		S1Top = 0.;
< 		S1Bot = 0.;
< 		// Check interaction type to deduce time constant
< 		if ( species == NR ) stRatio = 1.6;
< 		else stRatio = 0.6;
<
< 		for ( int i = 0; i < (int) floor(numPHETop); i++ ) {
< 			double phe1 = RandomGen::rndm()->rand_gauss(1., fdetector->get_sPEres());
< 			if ( phe1 < 0.0 || RandomGen::rndm()->rand_uniform() > fdetector->get_sPEeff() ) phe1 = 0.;
< 			double phe2 = 0.;
< 			if ( RandomGen::rndm()->rand_uniform() < fdetector->get_P_dphe() ) {
< 				phe2 = RandomGen::rndm()->rand_gauss(1., fdetector->get_sPEres());
< 				if ( phe2 < 0.0 || RandomGen::rndm()->rand_uniform() > fdetector->get_sPEeff() ) phe2 = 0.;
< 			}
<
< 			if ( RandomGen::rndm()->rand_uniform() < (stRatio/(1.+stRatio)) ) tau = 3.1;
< 			else tau = 24.;
<
< 			// spEff has units of phd, but before correcting from VUV v1.1 to v1.2
< 			if ( (phe1+phe2)*vuv1_1to1_2 > spEff_Run4 && RandomGen::rndm()->rand_exponential(tau) < coinWind ) {
< 				S1Top += (phe1+phe2);
< 			}
< 		}
<
< 		for ( int i = 0; i < (int) floor(numPHEBot); i++ ) {
< 			double phe1 = RandomGen::rndm()->rand_gauss(1., fdetector->get_sPEres());
< 			if ( phe1 < 0.0 || RandomGen::rndm()->rand_uniform() > fdetector->get_sPEeff() ) phe1 = 0.;
< 			double phe2 = 0.;
< 			if ( RandomGen::rndm()->rand_uniform() < fdetector->get_P_dphe() ) {
< 				phe2 = RandomGen::rndm()->rand_gauss(1., fdetector->get_sPEres());
< 				if ( phe2 < 0.0 || RandomGen::rndm()->rand_uniform() > fdetector->get_sPEeff() ) phe2 = 0.;
< 			}
<
< 			if ( RandomGen::rndm()->rand_uniform() < (stRatio/(1.+stRatio)) ) tau = 3.1;
< 			else tau = 24.;
<
< 			// spEff has units of phd, but before correcting from VUV v1.1 to v1.2
< 			if ( (phe1+phe2)*vuv1_1to1_2 > spEff_Run4 && RandomGen::rndm()->rand_exponential(tau) < coinWind ) {
< 				S1Bot += (phe1+phe2);
< 			}
< 		}
< 		// Correct for double-PE effect
< 		S1Top /= (1.+fdetector->get_P_dphe());
< 		S1Bot /= (1.+fdetector->get_P_dphe());
< 	}
< 	// If there is no coincidence requirement, simply calculate fluctuated number of PHE in top and bottom
< 	else {
< 		numPHETop += BinomFluct(numPHETop, fdetector->get_P_dphe());
< 		numPHEBot += BinomFluct(numPHEBot, fdetector->get_P_dphe());
< 		S1Top = RandomGen::rndm()->rand_gauss(numPHETop, fdetector->get_sPEres()*sqrt(numPHETop))/(1.+fdetector->get_P_dphe());
< 		S1Bot = RandomGen::rndm()->rand_gauss(numPHEBot, fdetector->get_sPEres()*sqrt(numPHEBot))/(1.+fdetector->get_P_dphe());
< 	}
<
< 	// Implement spike count correction based on lookup tables (TOP array)
< 	if (numSpikesTop < 81. && numSpikesTop > 0.5) {
< 		spikeCountCorrectionTop = SC_offset[numSpikesTop-1] + SC_lin[numSpikesTop-1]*S1Top +
< 										SC_quad[numSpikesTop-1]*pow(S1Top,2.) + SC_cubic[numSpikesTop-1]*pow(S1Top,3.);
< 	}
< 	else spikeCountCorrectionTop = 0.;
<
< 	// Implement spike count correction based on lookup tables (BOTTOM array)
< 	if (numSpikesBot < 81. && numSpikesBot > 0.5) {
< 		spikeCountCorrectionBot = SC_offset[numSpikesBot-1] + SC_lin[numSpikesBot-1]*S1Bot +
< 										SC_quad[numSpikesBot-1]*pow(S1Bot,2.) + SC_cubic[numSpikesBot-1]*pow(S1Bot,3.);
< 	}
< 	else spikeCountCorrectionBot = 0.;
<
< 	// Output variables: "cluged" corrected spike count, with and then without position corrections
< 	S1spike_OLc = (spikeCountCorrectionTop*double(numSpikesTop) + spikeCountCorrectionBot*double(numSpikesBot))*1.015;
< 	//newSpike[0] = (double) S1spike_raw; // use if you want the uncorrected spike count to NOT include Tomasz's correction
< 	newSpike[0] = S1spike_OLc; // non-position corrected (but spike-corrected) spike count
<   newSpike[1] = S1spike_OLc / fdetector->FitS1( dx, dy, dz ) * fdetector->FitS1( 0., 0., fdetector->get_TopDrift() - dS_mid * fdetector->get_dtCntr() ); // with position corrections
<
< 	return newSpike; // regular and position-corrected spike counts returned
< }
<

diff LibNEST/testNEST.cpp nest/testNEST.cpp
19,20d18
< #include "LUX_Run3.hh"
< #include "LUX_Run4.hh"
35,47c33,38
<   // Instantiate your own VDetector class here, then load into NEST class constructor
< 	//DetectorExample_XENON10* detector = new DetectorExample_XENON10();
< 	//LUX_Run3* detector = new LUX_Run3();
< 	LUX_Run4* detector = new LUX_Run4();
<
< 	// Run 3 (only uncomment if you need to change from default -> calibration)
< 	//detector->SetParameters_CH3T();
< 	//detector->SetParameters_DD();
< 	// Run 4 (will need to choose a time bin)
< 	detector->SetParameters_TimeBin1();
< 	//detector->SetParameters_TimeBin2();
< 	//detector->SetParameters_TimeBin3();
< 	//detector->SetParameters_TimeBin4();
---
>   // Instantiate your own VDetector class here, then load into NEST class
>   // constructor
>   DetectorExample_XENON10* detector = new DetectorExample_XENON10();
>
>   // Custom parameter modification functions
>   // detector->ExampleFunction();
```

