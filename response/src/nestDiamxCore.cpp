// LZ_model.cpp: main program to be used with python bindings
// This is a simplified version of execNEST.cpp
#include "analysis.hh"
#include "nestDiamxCore.hh"
#include "NEST.hh"
#include "LZ_WS2022.hh"
#include "LZ_WS2024.hh"

NEST::INTERACTION_TYPE getNESTType(string type) {
    if (type == "ER" || type == "beta" || type == "electron") {
        return NEST::INTERACTION_TYPE::beta;
    }
    else if (type == "gamma" || type == "photon" ||
             type == "DEC" || type == "EC") {
        return NEST::INTERACTION_TYPE::gammaRay;
    }
    else if (type == "NR") {
        return NEST::INTERACTION_TYPE::NR;
    }
    else {
        throw runtime_error("Unsupported type " + type);
    }
}

InferenceObservableArray LZModel(
    string detector_name,
    string type,
    uint64_t seed,
    string spectrumFileName,
    uint64_t numEvts,
    double dec_quenching_factor,
    double er_omega_scale,
    int detector_ablation_mode
) {
    if (er_omega_scale < 0.)
        throw runtime_error("er_omega_scale must be non-negative.");
    // Overriding default verbosity and maxS2
    int verbosity_lzmodel = -1;
    // Plot-only build for isotope peaks above the published WIMP-likelihood
    // S2 ceiling.  Keep the production DIAMX checkout unchanged.
    maxS1 = 1e7;
    maxS2 = 1e7;

    InferenceObservableArray InferenceObservable;

    // Detector initialization
    VDetector* detector = nullptr;
    vector<double> NRYieldsParam, ERYieldsParam, NRERWidthsParam;
    if (detector_name == "lz_ws2022" || detector_name == "lz_sr1") {
        detector = new LZ_Detector_2022();
        NRYieldsParam = default_NRYieldsParam;
        ERYieldsParam = default_ERYieldsParam;
        NRERWidthsParam = {0.4, 0.4, 0.04, 0.5, 0.19, 2.25, -0.0015, 0.046452, 0.205, 0.45, -0.2};
    }
    else if (detector_name == "lz_ws2024") {
        auto* lz_eft = new LZ_Detector_2024();
        lz_eft->SetEFT2024Configuration();
        detector = lz_eft;
        NRYieldsParam = lz_eft->get_nr_yield_params();
        ERYieldsParam = lz_eft->get_er_yield_params();
        NRERWidthsParam = lz_eft->get_nr_er_width_params();
    }
    else {
        throw runtime_error("Unsupported experiment " + detector_name);
    }
    // Diagnostic-only detector-response ablations. These do not alter the
    // production default (mode 0). Bit 0 removes the S1 single-PE/double-PE
    // broadening; bit 1 removes the S2 single-electron gain broadening.
    if (detector_ablation_mode & 1) {
        detector->set_sPEres(0.);
        detector->set_P_dphe(0.);
        detector->set_noiseLinear(0., detector->get_noiseLinear()[1]);
        detector->set_noiseQuadratic(0., detector->get_noiseQuadratic()[1]);
    }
    if (detector_ablation_mode & 2) {
        detector->set_s2Fano(0.);
        detector->set_noiseLinear(detector->get_noiseLinear()[0], 0.);
        detector->set_noiseQuadratic(detector->get_noiseQuadratic()[0], 0.);
    }
    // Scale only the ER non-binomial recombination-width amplitude.  All
    // yields, detector gains, binomial recombination, and detector smearing
    // remain at their selected LZ configuration values.
    NRERWidthsParam[7] *= er_omega_scale;
    if (NRERWidthsParam.size() > 13) NRERWidthsParam[13] *= er_omega_scale;

    // NEST calculation object
    NEST::NESTcalc n(detector);

    // Detector parameters
    double rho = n.SetDensity(detector->get_T_Kelvin(), detector->get_p_bar());
    double atomNum = 0; // Not used
    double massNum = detector->get_molarMass();
    // Calculate and print g1, g2 parameters (once per detector)
    vector<double> g2_params = n.CalculateG2(verbosity_lzmodel);
    double g2 = std::abs(g2_params[3]);
    double g1 = detector->get_g1();

    // Set seed
    if (seed == 0) {
        RandomGen::rndm()->SetSeed(time(nullptr));
    }
    else {
        RandomGen::rndm()->SetSeed(seed);
    }

    // Type translation
    NEST::INTERACTION_TYPE nestType = getNESTType(type);

    // Initialize energy sampler
    double defaultKeV = 0;
    EnergySampler* energySampler = nullptr;
    if (filesystem::exists(spectrumFileName)) {
        energySampler = new EnergySampler(spectrumFileName);
    }
    else {
        try {
            defaultKeV = stod(spectrumFileName);
        } catch (const std::exception& e) {
            throw runtime_error(
                "Invalid spectrumFileName. " + spectrumFileName + "can't be interpreted "
                + "as a valid file or a constant. \n" + "Error message: " + e.what()
            );
        }
        if (defaultKeV <= 0) {
            throw runtime_error("Invalid energy: " + to_string(defaultKeV) + "keV");
        }
    }

    // Event loop
    double keV = 0;
    double pos_x, pos_y, pos_z, r, phi;
    double field, vD, vD_middle;
    double driftTime;
    for (uint64_t j = 0; j < numEvts; ++j) {
        // Sample energy
        double keV = 0;
        if (energySampler) {
            keV = energySampler->sampleEnergy();
        }
        else {
            keV = defaultKeV;
        }
        if (keV <= 0) {
            // This should not happen, unless the spectrum is invalid
            throw runtime_error("Negative energy found. Please check the spectrum input.");
        }

        // Sample position
        bool z_new_flag = true; // I am never a fan of GOTO, so rewriting this part 
        while (z_new_flag) {
            pos_z = 0. + (detector->get_TopDrift() - 0.) *
                            RandomGen::rndm()->rand_uniform();  // initial guess
            r = detector->get_radius() * sqrt(RandomGen::rndm()->rand_uniform());
            phi = 2. * M_PI * RandomGen::rndm()->rand_uniform();
            pos_x = r * cos(phi);
            pos_y = r * sin(phi);

            // Electric field & drift velocity
            field = detector->FitEF(pos_x, pos_y, pos_z);
            vD = n.SetDriftVelocity(detector->get_T_Kelvin(), rho, field,
                                   detector->get_p_bar());
            vD_middle = vD;

            driftTime = (detector->get_TopDrift() - pos_z) / vD;
            // Reject positions outside the configured drift-time fiducial range.
            // The previous conjunction was inverted: it rejected valid positions
            // and accepted the first invalid one.
            z_new_flag = (driftTime < detector->get_dt_min()) ||
                         (driftTime > detector->get_dt_max());
        }
        
        // Yields & quanta generation
        NEST::YieldResult yields;
        NEST::QuantaResult quanta;
        NEST::NESTresult result;
        yields = n.GetYields(nestType, keV, rho, field, double(massNum),
                             double(atomNum), NRYieldsParam, ERYieldsParam);
        if (type == "DEC" || type == "EC") {
            double Nq = yields.ElectronYield + yields.PhotonYield;
            yields.ElectronYield *= dec_quenching_factor;
            yields.PhotonYield = Nq - yields.ElectronYield;
        }
        quanta = n.GetQuanta(yields, rho, NRERWidthsParam, 0.); // Turn off charge-distribution skewness only

        if (detector->get_noiseBaseline()[2] != 0. || detector->get_noiseBaseline()[3] != 0.)
            quanta.electrons += int(floor(
                RandomGen::rndm()->rand_gauss(
                    detector->get_noiseBaseline()[2],
                    detector->get_noiseBaseline()[3],
                    false
                ) + 0.5
            ));
        
        // Smeared position, not implemented but kept for extensibility
        double truthPos[3] = {pos_x, pos_y, pos_z};
        double smearPos[3] = {pos_x, pos_y, pos_z};
        double Nphd_S2 =
            g2 * quanta.electrons * exp(-driftTime / detector->get_eLife_us());
        if (!MCtruthPos && Nphd_S2 > PHE_MIN) {
            vector<double> xySmeared(2);
            xySmeared = n.xyResolution(pos_x, pos_y, Nphd_S2);
            smearPos[0] = xySmeared[0];
            smearPos[1] = xySmeared[1];
        }

        vector<int64_t> wf_time;
        vector<double> wf_amp;
        vector<double> scint =
          n.GetS1(quanta, pos_x, pos_y, pos_z, smearPos[0],
                  smearPos[1], smearPos[2], vD, vD_middle, nestType, j, field,
                  keV, s1CalculationMode, verbosity_lzmodel, wf_time, wf_amp);
        if (truthPos[2] < detector->get_cathode())
            quanta.electrons = 0;
        vector<double> scint2 = n.GetS2(quanta.electrons, truthPos[0], truthPos[1], truthPos[2],
                                        smearPos[0], smearPos[1], smearPos[2], driftTime, vD, j,
                                        field, s2CalculationMode, verbosity, wf_time, wf_amp,
                                        g2_params);
        
        if (std::abs(scint[5]) > minS1 && scint[5] < maxS1)
            InferenceObservable.s1c_phd.push_back(scint[5]);
        else
            InferenceObservable.s1c_phd.push_back(-999.);
        if (std::abs(scint2[7]) > minS2 && scint2[7] < maxS2)
            InferenceObservable.s2c_phd.push_back(scint2[7]);
        else
            InferenceObservable.s2c_phd.push_back(-999.);
        
        // Energy reconstruction
        double energy_rec = 0.0;
        double Nph = std::abs(scint[7]) / g1;
        double Ne = std::abs(scint2[7]) / g2;
        double Wq_eV = n.WorkFunction(rho, detector->get_molarMass(), detector->get_OldW13eV()).Wq_eV;
        if (yields.Lindhard > DBL_MIN && Nph > 0. && Ne > 0.) {
            if (ValidityTests::nearlyEqual(yields.Lindhard, 1.))
                energy_rec = (Nph + Ne) * Wq_eV * 1e-3;
            else {
                energy_rec = pow((Ne + Nph) / NRYieldsParam[0], 1. / NRYieldsParam[1]);
                Ne *= 1. - 1. / pow(1. + pow((energy_rec / NRYieldsParam[5]), NRYieldsParam[6]), NRYieldsParam[10]);
                Nph *= 1. - 1. / pow(1. + pow((energy_rec / NRYieldsParam[7]), NRYieldsParam[8]), NRYieldsParam[11]);
                energy_rec = pow((Ne + Nph) / NRYieldsParam[0], 1. / NRYieldsParam[1]);
            } 
        } else
            energy_rec = 0.;
        InferenceObservable.energy_rec.push_back(energy_rec);
    }

    return InferenceObservable;
}

InferenceObservableArray LZCascadeModel(
    string detector_name,
    vector<string> types,
    vector<double> energies,
    uint64_t seed,
    uint64_t numEvts,
    double er_omega_scale
) {
    if (er_omega_scale < 0.)
        throw runtime_error("er_omega_scale must be non-negative.");
    if (types.empty() || types.size() != energies.size())
        throw runtime_error("Cascade types and energies must have equal non-zero length.");
    for (double energy : energies)
        if (energy <= 0.) throw runtime_error("Cascade energies must be positive.");

    int verbosity_lzmodel = -1;
    // The cascade interface is used for the high-energy isotope test.  Match
    // LZModel's enlarged observable bounds instead of inheriting DIAMX's
    // low-energy maxS1=80 phd default.
    maxS1 = 1e7;
    maxS2 = 1e7;
    InferenceObservableArray output;

    VDetector* detector = nullptr;
    vector<double> NRYieldsParam, ERYieldsParam, NRERWidthsParam;
    if (detector_name == "lz_ws2022" || detector_name == "lz_sr1") {
        detector = new LZ_Detector_2022();
        NRYieldsParam = default_NRYieldsParam;
        ERYieldsParam = default_ERYieldsParam;
        NRERWidthsParam = {0.4, 0.4, 0.04, 0.5, 0.19, 2.25, -0.0015, 0.046452, 0.205, 0.45, -0.2};
    } else if (detector_name == "lz_ws2024") {
        auto* lz_eft = new LZ_Detector_2024();
        lz_eft->SetEFT2024Configuration();
        detector = lz_eft;
        NRYieldsParam = lz_eft->get_nr_yield_params();
        ERYieldsParam = lz_eft->get_er_yield_params();
        NRERWidthsParam = lz_eft->get_nr_er_width_params();
    } else {
        throw runtime_error("Unsupported experiment " + detector_name);
    }
    NRERWidthsParam[7] *= er_omega_scale;
    if (NRERWidthsParam.size() > 13) NRERWidthsParam[13] *= er_omega_scale;

    NEST::NESTcalc n(detector);
    double rho = n.SetDensity(detector->get_T_Kelvin(), detector->get_p_bar());
    double atomNum = 0.;
    double massNum = detector->get_molarMass();
    vector<double> g2_params = n.CalculateG2(verbosity_lzmodel);
    double g2 = std::abs(g2_params[3]);
    double g1 = detector->get_g1();
    RandomGen::rndm()->SetSeed(seed == 0 ? time(nullptr) : seed);

    double totalEnergy = 0.;
    for (double energy : energies) totalEnergy += energy;

    for (uint64_t j = 0; j < numEvts; ++j) {
        double pos_x, pos_y, pos_z, r, phi, field, vD, vD_middle, driftTime;
        bool retry = true;
        while (retry) {
            pos_z = detector->get_TopDrift() * RandomGen::rndm()->rand_uniform();
            r = detector->get_radius() * sqrt(RandomGen::rndm()->rand_uniform());
            phi = 2. * M_PI * RandomGen::rndm()->rand_uniform();
            pos_x = r * cos(phi); pos_y = r * sin(phi);
            field = detector->FitEF(pos_x, pos_y, pos_z);
            vD = n.SetDriftVelocity(detector->get_T_Kelvin(), rho, field,
                                   detector->get_p_bar());
            vD_middle = vD;
            driftTime = (detector->get_TopDrift() - pos_z) / vD;
            retry = (driftTime < detector->get_dt_min()) ||
                    (driftTime > detector->get_dt_max());
        }

        NEST::QuantaResult combined = {0, 0, 0, 0, 0., 0.};
        for (size_t k = 0; k < energies.size(); ++k) {
            NEST::INTERACTION_TYPE interaction = getNESTType(types[k]);
            NEST::YieldResult yields = n.GetYields(
                interaction, energies[k], rho, field, massNum, atomNum,
                NRYieldsParam, ERYieldsParam);
            NEST::QuantaResult part = n.GetQuanta(
                yields, rho, NRERWidthsParam, 0.);
            combined.photons += part.photons;
            combined.electrons += part.electrons;
            combined.ions += part.ions;
            combined.excitons += part.excitons;
        }

        if (detector->get_noiseBaseline()[2] != 0. || detector->get_noiseBaseline()[3] != 0.)
            combined.electrons += int(floor(RandomGen::rndm()->rand_gauss(
                detector->get_noiseBaseline()[2], detector->get_noiseBaseline()[3], false) + 0.5));

        double truthPos[3] = {pos_x, pos_y, pos_z};
        double smearPos[3] = {pos_x, pos_y, pos_z};
        double Nphd_S2 = g2 * combined.electrons * exp(-driftTime / detector->get_eLife_us());
        if (!MCtruthPos && Nphd_S2 > PHE_MIN) {
            vector<double> xySmeared = n.xyResolution(pos_x, pos_y, Nphd_S2);
            smearPos[0] = xySmeared[0]; smearPos[1] = xySmeared[1];
        }

        vector<int64_t> wf_time;
        vector<double> wf_amp;
        NEST::INTERACTION_TYPE displayType = getNESTType(types[0]);
        vector<double> scint = n.GetS1(
            combined, pos_x, pos_y, pos_z, smearPos[0], smearPos[1], smearPos[2],
            vD, vD_middle, displayType, j, field, totalEnergy,
            s1CalculationMode, verbosity_lzmodel, wf_time, wf_amp);
        if (truthPos[2] < detector->get_cathode()) combined.electrons = 0;
        vector<double> scint2 = n.GetS2(
            combined.electrons, truthPos[0], truthPos[1], truthPos[2],
            smearPos[0], smearPos[1], smearPos[2], driftTime, vD, j, field,
            s2CalculationMode, verbosity, wf_time, wf_amp, g2_params);

        output.s1c_phd.push_back(
            std::abs(scint[5]) > minS1 && scint[5] < maxS1 ? scint[5] : -999.);
        output.s2c_phd.push_back(
            std::abs(scint2[7]) > minS2 && scint2[7] < maxS2 ? scint2[7] : -999.);
        double Nph = std::abs(scint[7]) / g1;
        double Ne = std::abs(scint2[7]) / g2;
        double Wq_eV = n.WorkFunction(rho, detector->get_molarMass(), detector->get_OldW13eV()).Wq_eV;
        output.energy_rec.push_back((Nph > 0. && Ne > 0.) ? (Nph + Ne) * Wq_eV * 1e-3 : 0.);
    }
    delete detector;
    return output;
}
