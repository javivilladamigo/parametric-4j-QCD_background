#include <TMath.h>
#include <TRandom2.h>
#include <TF1.h>
#include <TLorentzVector.h>
#include <TH1D.h>
#include <TCanvas.h>

struct jetEvent_t {
    Double_t pt;
    Double_t eta;
    Double_t phi;
};

/////////////////////////////////////////////////////////////////////////////////////////////////
// \author: Javier Mariño Villadamigo                                                          //
//                                                                                             //
// Date: 30/11/2022                                                                            //
//                                                                                             //
// Thesis: parametric model for the generation of a 4-jets QCD background                      //
//
//
//
//
//
// NOTE:
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////


/*
Double_t velocity(Double_t Energy, Double_t x1, Double_t x2, Double_t m) {
    return Energy * (x1 - x2) / TMath::Sqrt(Energy * Energy * TMath::Power(x1-x2, 2) + 4*m*m);
}
*/

Double_t gamma_v(Double_t v) {
    // gamma factor
    return 1 / TMath::Sqrt(1 - v * v);
}

Double_t transversal_coordinates_tfm(Double_t p_x_or_y_CoM, Double_t pz_CoM, Double_t v, Double_t m) {
    // transformation of the transverse components of 4-momentum (in CoM) to the LAB moving with velocity -v wrt to CoM along z
    Double_t gamma_factor = gamma_v(v);
    Double_t p_x_or_y_LAB = p_x_or_y_CoM / (gamma_factor * (gamma_factor * m + v * pz_CoM));
    return p_x_or_y_LAB;
}

Double_t longitudinal_coordinate_tfm(Double_t pz_CoM, Double_t v, Double_t m) {
    // transformation of the longitudinal components of 4-momentum
    Double_t gamma_factor = gamma_v(v);
    Double_t pz_LAB = (pz_CoM + gamma_factor * m * v) / (gamma_factor * m + v * pz_CoM);
    return pz_LAB;
}

Double_t energy_tfm(Double_t pz_CoM, Double_t v, Double_t E_CoM) {
    // relativistic energy shift
    Double_t gamma_factor = gamma_v(v);
    Double_t E_LAB = gamma_factor * ( E_CoM - v * pz_CoM);
    return E_LAB;
}

TLorentzVector boost_toLAB(TLorentzVector p_CoM, Double_t v_CoM, Double_t m) {
    // Lorentz transformation of 4-momenta from CoM to LAB frame
    TLorentzVector p_LAB;
    p_LAB.SetPxPyPzE(transversal_coordinates_tfm(p_CoM.Px(), p_CoM.Pz(), v_CoM, m), transversal_coordinates_tfm(p_CoM.Py(), p_CoM.Pz(), v_CoM, m), longitudinal_coordinate_tfm(p_CoM.Pz(), v_CoM, m), energy_tfm(p_CoM.Pz(), v_CoM, p_CoM.E()));
    return p_LAB; 
}

Double_t thetaLAB_from_pLAB(TLorentzVector p_LAB) {
    // Lorentz transformation of the angle of the particle motion relative to the z axis in a coordinate frame
    // arccos (pz / |p|)
    return TMath::ACos(p_LAB.Pz() / TMath::Sqrt(TMath::Power(p_LAB.Px(), 2) + TMath::Power(p_LAB.Py(), 2) + TMath::Power(p_LAB.Pz(), 2)));
}

Double_t eta_from_theta(Double_t theta) {
    return (-TMath::Log(TMath::Tan(theta / 2)));
}



void parametric_4j_QCD_background(Int_t Nevents = 10000) {

    delete gRandom;
    //int seed = 0;
    gRandom = new TRandom2();

    

    //TNtuple *mytuple = new TNtuple("mytuple", "parametric_4j_QCD_background", "j_Px:j_Py:j_Pz:j_E:j_eta:j_phi");
    //TFile *outputFile = new TFile("NTuple.root", "recreate");

    // VARIABLES
        // pre-jets
        Double_t tau = 1./10.;                                                                              // Bjorken scaling slope
        Double_t E_beam = 6.5;                                                                              // TeV of energy for each proton
        Double_t mass = 2E-6;                                                                               // 2 MeV of mass assumed for the quark (in TeV)
        Double_t E_i1_CoM, E_i2_CoM;                                                                        // TeV of energy for the incoming partons (in CoM frame)
        Double_t sqrt_s;                                                                                    // √s 
        Double_t E1_prime_CoM, E2_prime_CoM;                                                                // TeV of energy for the outgoing partons (in CoM frame)
        Double_t E1_prime_LAB, E2_prime_LAB;                                                                // TeV of energy for the outgoing partons (in LAB frame)
        Double_t theta_CoM, phi1, phi2, theta1_LAB, theta2_LAB;                                             // angles of the outgoing partons
        TLorentzVector p1_prime_CoM, p2_prime_CoM, pSUM_CoM, p1_prime_LAB, p2_prime_LAB;                    // Lorentz vectors of partons after collision
        Double_t velocity_CoM_frame;                                                                        // velocity of CoM wrt to LAB
        
    
        // post-jets
        Double_t x1, x2, split_prob1, split_prob2;                                                          // Bjorken scalings and probability of jet splitting
        Double_t E_j1_CoM, E_j2_CoM, E_j3_CoM, E_j4_CoM;                                                    // TeV of energy of the jets (in CoM frame)
        Double_t E_j1_LAB, E_j2_LAB, E_j3_LAB, E_j4_LAB;                                                    // TeV of energy of the jets (in CoM frame)
        Double_t theta_j1_CoM, theta_j2_CoM, theta_j3_CoM, theta_j4_CoM, phi_j1_CoM, phi_j2_CoM, phi_j3_CoM, phi_j4_CoM;    // angles of the jets emitted in the CoM; j2 and j4 are collinear with j1 and j3 respectively
        Double_t theta_j1_LAB, theta_j2_LAB, theta_j3_LAB, theta_j4_LAB, phi_j1_LAB, phi_j2_LAB, phi_j3_LAB, phi_j4_LAB;    // angles of the jets emitted in the CoM; j2 and j4 are collinear with j1 and j3 respectively
        TLorentzVector p_j1_CoM, p_j2_CoM, p_j3_CoM, p_j4_CoM;                                              // Lorentz vectors of jets (in CoM frame)
        TLorentzVector p_j1_LAB, p_j2_LAB, p_j3_LAB, p_j4_LAB;                                              // Lorentz vectors of jets (in LAB frame)
        
        
        TF1 *F = new TF1("F", "0.5 + 2 * (TMath::Erf(50-x) + TMath::Erf(30))");                             // detection efficiency
        TF1 *B = new TF1("B", "0.7 * (0.5 + 2 * TMath::Erf((1.8 - TMath::Abs(x)) / 0.3))");                 // b-tagging efficiency
        Double_t F_j1, F_j2, F_j3, F_j4;                                                                    // uniformly random numbers to be compared to F distr (detection efficiency)
        Double_t F_th;                                                                                      // upper random value for all the F_j's
        Double_t B_j1, B_j2, B_j3, B_j4;                                                                    // uniformly random numbers to be compared to B distr (b-tagging efficiency)
        Double_t B_th;                                                                                      // upper random value for all the B_j's
        Int_t Nbtaggedjets;                                                                                 // nb of jets b-tagged


    // HISTOGRAMS
        TH1D *hist_x1 = new TH1D("", "Bjorken x1; x1; Events", 200, 0, 2);
        TH1D *hist_x2 = new TH1D("", "Bjorken x2; x2; Events", 200, 0, 2);
        

        TH1D *hist_E1_prime_CoM = new TH1D("", "E parton 1 (blue CoM and red LAB); TeV; Events", 250, -1, 20);
        hist_E1_prime_CoM->SetLineColor(kBlue);
        TH1D *hist_E1_prime_LAB = new TH1D("", "", 250, -1, 20);
        hist_E1_prime_LAB->SetLineColor(kRed);


        TH1D *hist_E2_prime_CoM = new TH1D("", "E parton 2 (blue CoM and red LAB); TeV; Events", 250, -1, 20);
        hist_E2_prime_CoM->SetLineColor(kBlue);
        TH1D *hist_E2_prime_LAB = new TH1D("", "", 250, -1, 20);
        hist_E2_prime_LAB->SetLineColor(kRed);


        TH1D *hist_sumE_CoM = new TH1D("", "sumE of 2partons (blue CoM and red LAB); TeV; Events", 250, -1, 40);
        hist_sumE_CoM->SetLineColor(kBlue);
        TH1D *hist_sumE_LAB = new TH1D("", "", 250, -1, 40);
        hist_sumE_LAB->SetLineColor(kRed);


        TH1D *hist_v = new TH1D("", "Velocity of the CoM frame (in c units); v/c; Events", 250, -1.2, 1.2);
        hist_v->SetLineColor(kBlue);


        TH1D *hist_thetaLAB1 = new TH1D("", "thetaLAB (blue parton1 and red parton2); #theta(rad); Events", 250, -2, TMath::Pi() + 1);
        hist_thetaLAB1->SetLineColor(kBlue);
        TH1D *hist_thetaLAB2 = new TH1D("", "", 250, -2, TMath::Pi() + 1);
        hist_thetaLAB2->SetLineColor(kRed);


        TH1D *hist_E4j_CoM = new TH1D("", "energy sum of 4jets (blue CoM and red LAB); TeV; Events", 250, -1., 40.);
        hist_E4j_CoM->SetLineColor(kBlue);
        TH1D *hist_E4j_LAB = new TH1D("", "", 250, -1., 40.);
        hist_E4j_LAB->SetLineColor(kRed);


        TH1D *hist_thetaj1_LAB = new TH1D("", "theta LAB 4jets; #theta(rad); Events", 250, -10., 10.);
        hist_thetaj1_LAB->SetLineColor(kBlue);
        TH1D *hist_thetaj2_LAB = new TH1D("", "", 250, -10., 10.);
        hist_thetaj2_LAB->SetLineColor(kRed);
        TH1D *hist_thetaj3_LAB = new TH1D("", "", 250, -10., 10.);
        hist_thetaj3_LAB->SetLineColor(kGreen);
        TH1D *hist_thetaj4_LAB = new TH1D("", "", 250, -10., 10.);
        hist_thetaj4_LAB->SetLineColor(kBlack);
        

        TH1D *hist_phij1_LAB = new TH1D("", "phi LAB 4jets; #phi(rad); Events", 250, -10., 10.);
        hist_phij1_LAB->SetLineColor(kBlue);
        TH1D *hist_phij2_LAB = new TH1D("", "", 250, -10., 10.);
        hist_phij2_LAB->SetLineColor(kRed);
        TH1D *hist_phij3_LAB = new TH1D("", "", 250, -10., 10.);
        hist_phij3_LAB->SetLineColor(kGreen);
        TH1D *hist_phij4_LAB = new TH1D("", "", 250, -10., 10.);
        hist_phij4_LAB->SetLineColor(kBlack);


        TH1D *hist_etaj1 = new TH1D("", "eta 4jets; #eta; Events", 250, -10., 10.);
        hist_etaj1->SetLineColor(kBlue);
        TH1D *hist_etaj2 = new TH1D("", "", 250, -10., 10.);
        hist_etaj2->SetLineColor(kRed);
        TH1D *hist_etaj3 = new TH1D("", "", 250, -10., 10.);
        hist_etaj3->SetLineColor(kGreen);
        TH1D *hist_etaj4 = new TH1D("", "", 250, -10., 10.);
        hist_etaj4->SetLineColor(kBlack);


        TH1D *hist_pTj1 = new TH1D("", "pT 4jets; p_{T} (TeV); Events", 250, -0.1, 30);
        hist_pTj1->SetLineColor(kBlue);
        TH1D *hist_pTj2 = new TH1D("", "", 250, -0.1, 30);
        hist_pTj2->SetLineColor(kRed);
        TH1D *hist_pTj3 = new TH1D("", "", 250, -0.1, 30);
        hist_pTj3->SetLineColor(kGreen);
        TH1D *hist_pTj4 = new TH1D("", "", 250, -0.1, 30);
        hist_pTj4->SetLineColor(kBlack);


        TH1D *hist_F = new TH1D("", "F", 250, 0., 1.);
        TH1D *hist_B = new TH1D("", "B", 250, 0., 1.);

    
    // create output file and output TTree
    TString outfname;
    outfname.Form("4jQCDbackground_%ievents.root", Nevents);
    TFile *outFile = new TFile(outfname, "recreate");
    TTree *Events_4jTree = new TTree("Events_4j", "Events_4j");
    
    jetEvent_t jetEvent1, jetEvent2, jetEvent3, jetEvent4;
    TBranch *j1 = Events_4jTree->Branch("j1", &jetEvent1);
    TBranch *j2 = Events_4jTree->Branch("j2", &jetEvent2);
    TBranch *j3 = Events_4jTree->Branch("j3", &jetEvent3);
    TBranch *j4 = Events_4jTree->Branch("j4", &jetEvent4);




    // event generation
    for (Int_t i = 0; i < Nevents; i++)
    {   
        x1 = gRandom->Exp(tau); // exp (-t / tau)
        x2 = gRandom->Exp(tau); // exp (-t / tau)
        while (x1 > 1. || x2 > 1.)
        {
            x1 = gRandom->Exp(tau);
            x2 = gRandom->Exp(tau);
        }

        // calculate energy of the incoming partons
        E_i1_CoM = x1 * E_beam;
        E_i2_CoM = x2 * E_beam;
        
        // energy in the CoM
        sqrt_s = E_i1_CoM + E_i2_CoM; 

        // equal energy partition between partons in the CoM
        E1_prime_CoM = sqrt_s / 2.; 
        E2_prime_CoM = sqrt_s / 2.;

        // angles of the outgoing partons in the 2->2 process
        theta_CoM = gRandom->Uniform(0, TMath::Pi());
        phi1 = gRandom->Uniform(0, 2 * TMath::Pi()); 
        phi2 = phi1 + TMath::Pi();

        // considering 0 mass also for the outgoing partons
        p1_prime_CoM.SetPxPyPzE(E1_prime_CoM * TMath::Sin(theta_CoM) * TMath::Cos(phi1), E1_prime_CoM * TMath::Sin(theta_CoM) * TMath::Sin(phi1), E1_prime_CoM * TMath::Cos(theta_CoM), E1_prime_CoM);
        p2_prime_CoM.SetPxPyPzE(E2_prime_CoM * TMath::Sin(TMath::Pi() - theta_CoM) * TMath::Cos(phi2), E2_prime_CoM * TMath::Sin(TMath::Pi() - theta_CoM) * TMath::Sin(phi2), E2_prime_CoM * TMath::Cos(TMath::Pi() - theta_CoM), E2_prime_CoM); // theta for second parton is (pi - theta) in CoM
        
        // 4-momentum of the composite system in the CoM - 3-momentum should be ~0
        pSUM_CoM.SetPxPyPzE(p1_prime_CoM.Px() + p2_prime_CoM.Px(), p1_prime_CoM.Py() + p2_prime_CoM.Py(), p1_prime_CoM.Pz() + p2_prime_CoM.Pz(), p1_prime_CoM.E() + p2_prime_CoM.E());
        // Indeed the CoM of the composite components are nearly 0 except rounding errors
        // cout << p1_prime_CoM.Px() + p2_prime_CoM.Px() << " " << p1_prime_CoM.Py() + p2_prime_CoM.Py() << " " << p1_prime_CoM.Pz() + p2_prime_CoM.Pz() << " " << p1_prime_CoM.E() + p2_prime_CoM.E() << endl;
        
        // since I have now the energy in the CoM, and (x1-x2) will indicate the direction in which the CoM frame will move wrt the LAB
        // I can calculate the velocity of such CoM frame
        velocity_CoM_frame = x1 - x2; // c units

        // boosting momenta to the LAB frame
        p1_prime_LAB = boost_toLAB(p1_prime_CoM, velocity_CoM_frame, mass);
        p2_prime_LAB = boost_toLAB(p2_prime_CoM, velocity_CoM_frame, mass);

        // calculating energies in the LAB frame (just for visualization, they are implicitly calculated in the boost_toLAB)
        E1_prime_LAB = energy_tfm(p1_prime_CoM.Pz(), velocity_CoM_frame, p1_prime_CoM.E());
        E2_prime_LAB = energy_tfm(p2_prime_CoM.Pz(), velocity_CoM_frame, p2_prime_CoM.E());


        theta1_LAB = TMath::ACos(p1_prime_LAB.Pz() / TMath::Sqrt(TMath::Power(p1_prime_LAB.Px(), 2) + TMath::Power(p1_prime_LAB.Py(), 2) + TMath::Power(p1_prime_LAB.Pz(), 2)));
        theta2_LAB = TMath::ACos(p2_prime_LAB.Pz() / TMath::Sqrt(TMath::Power(p2_prime_LAB.Px(), 2) + TMath::Power(p2_prime_LAB.Py(), 2) + TMath::Power(p2_prime_LAB.Pz(), 2)));





        // histograms up to before jets generation
        hist_x1->Fill(x1);
        hist_x2->Fill(x2);

        hist_E1_prime_CoM->Fill(E1_prime_CoM);
        hist_E1_prime_LAB->Fill(E1_prime_LAB);
        hist_E2_prime_CoM->Fill(E2_prime_CoM);
        hist_E2_prime_LAB->Fill(E2_prime_LAB);

        hist_sumE_CoM->Fill(E1_prime_CoM + E2_prime_CoM);
        hist_sumE_LAB->Fill(E1_prime_LAB + E2_prime_LAB);

        hist_v->Fill(velocity_CoM_frame);

        hist_thetaLAB1->Fill(theta1_LAB);
        hist_thetaLAB2->Fill(theta2_LAB);






        // probabilities of splitting
        split_prob1 = gRandom->Uniform(-1, 1);
        split_prob2 = gRandom->Uniform(-1, 1);

        if (split_prob1 > 0 && split_prob2 > 0) // splitting to 4 partons
        {

            // will use split prob also as a separator of energy since in the case of splitting split_prob is in (0, 1)
            E_j1_CoM = split_prob1 * E1_prime_CoM; // using the energy in CoM !!
            E_j2_CoM = (1 - split_prob1) * E1_prime_CoM;
            E_j3_CoM = split_prob2 * E2_prime_CoM;
            E_j4_CoM = (1 - split_prob2) * E2_prime_CoM;

            hist_E4j_CoM->Fill(E_j1_CoM + E_j2_CoM + E_j3_CoM + E_j4_CoM);

            // theta angles of jets in the CoM
            theta_j1_CoM = gRandom->Uniform(0, TMath::Pi());
            theta_j2_CoM = TMath::Pi() - theta_j1_CoM;
            theta_j3_CoM = gRandom->Uniform(0, TMath::Pi());
            theta_j4_CoM = TMath::Pi() - theta_j3_CoM;

            // phi angles of jets in the CoM
            phi_j1_CoM = gRandom->Uniform(0, 2 * TMath::Pi());
            phi_j2_CoM = phi_j1_CoM + TMath::Pi();
            phi_j3_CoM = gRandom->Uniform(0, 2 * TMath::Pi());
            phi_j4_CoM = phi_j3_CoM + TMath::Pi();

            // CoM 4-momenta of jets
            p_j1_CoM.SetPxPyPzE(E_j1_CoM * TMath::Sin(theta_j1_CoM) * TMath::Cos(phi_j1_CoM), E_j1_CoM * TMath::Sin(theta_j1_CoM) * TMath::Sin(phi_j1_CoM), E_j1_CoM * TMath::Cos(theta_j1_CoM), E_j1_CoM);
            p_j2_CoM.SetPxPyPzE(E_j2_CoM * TMath::Sin(theta_j2_CoM) * TMath::Cos(phi_j2_CoM), E_j2_CoM * TMath::Sin(theta_j2_CoM) * TMath::Sin(phi_j2_CoM), E_j2_CoM * TMath::Cos(theta_j2_CoM), E_j2_CoM);
            p_j3_CoM.SetPxPyPzE(E_j3_CoM * TMath::Sin(theta_j3_CoM) * TMath::Cos(phi_j3_CoM), E_j3_CoM * TMath::Sin(theta_j3_CoM) * TMath::Sin(phi_j3_CoM), E_j3_CoM * TMath::Cos(theta_j3_CoM), E_j3_CoM);
            p_j4_CoM.SetPxPyPzE(E_j4_CoM * TMath::Sin(theta_j4_CoM) * TMath::Cos(phi_j4_CoM), E_j4_CoM * TMath::Sin(theta_j4_CoM) * TMath::Sin(phi_j4_CoM), E_j4_CoM * TMath::Cos(theta_j4_CoM), E_j4_CoM);

            
            // boosting jets 4-momenta to the LAB frame
            p_j1_LAB = boost_toLAB(p_j1_CoM, velocity_CoM_frame, mass / 2.);
            p_j2_LAB = boost_toLAB(p_j2_CoM, velocity_CoM_frame, mass / 2.);
            p_j3_LAB = boost_toLAB(p_j3_CoM, velocity_CoM_frame, mass / 2.);
            p_j4_LAB = boost_toLAB(p_j4_CoM, velocity_CoM_frame, mass / 2.);
            
            // thetaLAB of jets
            theta_j1_LAB = thetaLAB_from_pLAB(p_j1_LAB);
            theta_j2_LAB = thetaLAB_from_pLAB(p_j2_LAB);
            theta_j3_LAB = thetaLAB_from_pLAB(p_j3_LAB);
            theta_j4_LAB = thetaLAB_from_pLAB(p_j4_LAB);

            hist_thetaj1_LAB->Fill(theta_j1_LAB);
            hist_thetaj2_LAB->Fill(theta_j2_LAB);
            hist_thetaj3_LAB->Fill(theta_j3_LAB);
            hist_thetaj4_LAB->Fill(theta_j4_LAB);


            // pT and eta
            Double_t pT_j1, pT_j2, pT_j3, pT_j4;
            Double_t eta_j1, eta_j2, eta_j3, eta_j4;

            pT_j1 = TMath::Sqrt(p_j1_LAB.Px()*p_j1_LAB.Px() + p_j1_LAB.Py()*p_j1_LAB.Py()+p_j1_LAB.Pz()*p_j1_LAB.Pz()) * TMath::Sin(theta_j1_LAB);
            pT_j2 = TMath::Sqrt(p_j2_LAB.Px()*p_j2_LAB.Px() + p_j2_LAB.Py()*p_j2_LAB.Py()+p_j2_LAB.Pz()*p_j2_LAB.Pz()) * TMath::Sin(theta_j2_LAB);
            pT_j3 = TMath::Sqrt(p_j3_LAB.Px()*p_j3_LAB.Px() + p_j3_LAB.Py()*p_j3_LAB.Py()+p_j3_LAB.Pz()*p_j3_LAB.Pz()) * TMath::Sin(theta_j3_LAB);
            pT_j4 = TMath::Sqrt(p_j4_LAB.Px()*p_j4_LAB.Px() + p_j4_LAB.Py()*p_j4_LAB.Py()+p_j4_LAB.Pz()*p_j4_LAB.Pz()) * TMath::Sin(theta_j4_LAB);

            hist_pTj1->Fill(pT_j1);
            hist_pTj2->Fill(pT_j2);
            hist_pTj3->Fill(pT_j3);
            hist_pTj4->Fill(pT_j4);

            eta_j1 = - TMath::Log(TMath::Tan(theta_j1_LAB / 2));
            eta_j2 = - TMath::Log(TMath::Tan(theta_j2_LAB / 2));
            eta_j3 = - TMath::Log(TMath::Tan(theta_j3_LAB / 2));
            eta_j4 = - TMath::Log(TMath::Tan(theta_j4_LAB / 2));

            hist_etaj1->Fill(eta_j1);
            hist_etaj2->Fill(eta_j2);
            hist_etaj3->Fill(eta_j3);
            hist_etaj4->Fill(eta_j4);

            phi_j1_LAB = 2 * TMath::ATan(p_j1_LAB.Py() / p_j1_LAB.Px());
            phi_j2_LAB = 2 * TMath::ATan(p_j2_LAB.Py() / p_j2_LAB.Px());
            phi_j3_LAB = 2 * TMath::ATan(p_j3_LAB.Py() / p_j3_LAB.Px());
            phi_j4_LAB = 2 * TMath::ATan(p_j4_LAB.Py() / p_j4_LAB.Px());

            hist_phij1_LAB->Fill(phi_j1_LAB);
            hist_phij2_LAB->Fill(phi_j2_LAB);
            hist_phij3_LAB->Fill(phi_j3_LAB);
            hist_phij4_LAB->Fill(phi_j4_LAB);

            E_j1_LAB = energy_tfm(p_j1_CoM.Pz(), velocity_CoM_frame, p_j1_CoM.E());
            E_j2_LAB = energy_tfm(p_j2_CoM.Pz(), velocity_CoM_frame, p_j2_CoM.E());
            E_j3_LAB = energy_tfm(p_j3_CoM.Pz(), velocity_CoM_frame, p_j3_CoM.E());
            E_j4_LAB = energy_tfm(p_j4_CoM.Pz(), velocity_CoM_frame, p_j4_CoM.E());

            hist_E4j_LAB->Fill(E_j1_LAB + E_j2_LAB + E_j3_LAB + E_j4_LAB);

            p_j1_LAB.SetPtEtaPhiE(pT_j1, eta_j1, phi_j1_LAB, E_j1_LAB);
            p_j2_LAB.SetPtEtaPhiE(pT_j2, eta_j2, phi_j2_LAB, E_j2_LAB);
            p_j3_LAB.SetPtEtaPhiE(pT_j3, eta_j3, phi_j3_LAB, E_j3_LAB);
            p_j4_LAB.SetPtEtaPhiE(pT_j4, eta_j4, phi_j4_LAB, E_j4_LAB);



            // detection acceptance
            F_j1 = gRandom->Uniform(0, 1);
            F_j2 = gRandom->Uniform(0, 1);
            F_j3 = gRandom->Uniform(0, 1);
            F_j4 = gRandom->Uniform(0, 1);
            F_th = F->GetRandom();
            hist_F->Fill(F_th);

            if (F_j1 < F_th && F_j2 < F_th && F_j3 < F_th && F_j4 < F_th) // all 4 jets detected
            {
                // b-tagging efficiency
                B_j1 = gRandom->Uniform(0, 1);
                B_j2 = gRandom->Uniform(0, 1);
                B_j3 = gRandom->Uniform(0, 1);
                B_j4 = gRandom->Uniform(0, 1);
                B_th = F->GetRandom();

                if (B_j1 < B_th && B_j2 < B_th && B_j3 < B_th && B_j4 < B_th) // all 4 jets b-tagged
                {   
                    Nbtaggedjets = 4;
                    jetEvent1.pt = p_j1_LAB.Pt();
                    jetEvent1.eta = p_j1_LAB.Eta();
                    jetEvent1.phi = p_j1_LAB.Phi();

                    jetEvent2.pt = p_j2_LAB.Pt();
                    jetEvent2.eta = p_j2_LAB.Eta();
                    jetEvent2.phi = p_j2_LAB.Phi();

                    jetEvent3.pt = p_j3_LAB.Pt();
                    jetEvent3.eta = p_j3_LAB.Eta();
                    jetEvent3.phi = p_j3_LAB.Phi();

                    jetEvent4.pt = p_j4_LAB.Pt();
                    jetEvent4.eta = p_j4_LAB.Eta();
                    jetEvent4.phi = p_j4_LAB.Phi();
                    
                    
                    Events_4jTree->Fill();
                }

                if ((B_j1 < B_th && B_j2 < B_th && B_j3 < B_th && B_j4 > B_th) || (B_j1 < B_th && B_j2 < B_th && B_j3 > B_th && B_j4 < B_th) || (B_j1 < B_th && B_j2 > B_th && B_j3 < B_th && B_j4 < B_th) || ((B_j1 > B_th && B_j2 < B_th && B_j3 < B_th && B_j4 < B_th))) // 3 jets b-tagged
                {
                    Nbtaggedjets = 3;
                    //mytuple->Fill(p_j1_CoM.Px(), p_j1_CoM.Py(), p_j1_CoM.Pz(), p_j1_CoM.E(), p_j1_CoM.Eta(), p_j1_CoM.Phi(), p_j2_CoM.Px(), p_j2_CoM.Py(), p_j2_CoM.Pz(), p_j2_CoM.E(), p_j2_CoM.Eta(), p_j2_CoM.Phi(), p_j3_CoM.Px(), p_j3_CoM.Py(), p_j3_CoM.Pz(), p_j3_CoM.E(), p_j3_CoM.Eta(), p_j3_CoM.Phi(), p_j4_CoM.Px(), p_j4_CoM.Py(), p_j4_CoM.Pz(), p_j4_CoM.E(), p_j4_CoM.Eta(), p_j4_CoM.Phi(), Nbtaggedjets)
                }

            }
            else
            {
                continue;
            }

        }
        else if ((split_prob1 > 0 && split_prob2 < 0) || (split_prob1 < 0 && split_prob2 > 0)) // splitting to three partons
        {
            //cout << "three jets" << endl;
        }
        else if (split_prob1 < 0 && split_prob2 < 0) // no splitting (two partons)
        {
            //cout << "two jets" << endl;
        }
    }
    outFile->cd();
    Events_4jTree->Write("", TObject::kOverwrite);
    outFile->Close();
    delete outFile;

    // plotting for checking
    TCanvas *c1 = new TCanvas("c1", "x1");
    c1->SetLogy();
    hist_x1->Draw();

    TCanvas *c2 = new TCanvas("c2", "x2");
    c2->SetLogy();
    hist_x2->Draw();


    TCanvas *c6 = new TCanvas("c6", "Vel. CoM");
    c6->SetLogy();
    hist_v->Draw();


    TCanvas *c3 = new TCanvas("c3", "E1");
    c3->SetLogy();
    hist_E1_prime_CoM->Draw();
    hist_E1_prime_LAB->Draw("same");

    TCanvas *c4 = new TCanvas("c4", "E2");
    c4->SetLogy();
    hist_E2_prime_CoM->Draw();
    hist_E2_prime_LAB->Draw("same");


    TCanvas *c7 = new TCanvas("c7", "theta_LAB");
    hist_thetaLAB1->Draw();
    hist_thetaLAB2->Draw("same");


    TCanvas *c8 = new TCanvas("c8", "sumE 2partons");
    c8->SetLogy();
    hist_sumE_CoM->Draw();
    hist_sumE_LAB->Draw("same");


    TCanvas *c9 = new TCanvas("c9", "sumE 4jets");
    c9->SetLogy();
    hist_E4j_CoM->Draw("");
    hist_E4j_LAB->Draw("same");


    TCanvas *c10 = new TCanvas("c10", "thetaLAB 4jets");
    hist_thetaj1_LAB->Draw("");
    hist_thetaj2_LAB->Draw("same");
    hist_thetaj3_LAB->Draw("same");
    hist_thetaj4_LAB->Draw("same");

    TCanvas *c11 = new TCanvas("c11", "phiLAB 4jets");
    hist_phij1_LAB->Draw("");
    hist_phij2_LAB->Draw("same");
    hist_phij3_LAB->Draw("same");
    hist_phij4_LAB->Draw("same");

    TCanvas *c12 = new TCanvas("c12", "eta 4jets");
    hist_etaj1->Draw("");
    hist_etaj2->Draw("same");
    hist_etaj3->Draw("same");
    hist_etaj4->Draw("same");


    TCanvas *c13 = new TCanvas("c13", "pT 4jets");
    c12->SetLogy();
    hist_pTj1->Draw("");
    hist_pTj2->Draw("same");
    hist_pTj3->Draw("same");
    hist_pTj4->Draw("same");
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// using the Boost function from TLorentzVector
void parametric_4j_QCD_background_autoBoost(Int_t Nevents = 500000) {

    delete gRandom;
    //int seed = 0;
    gRandom = new TRandom2();

    

    // VARIABLES
        // pre-jets
        Double_t tau = 1./10.;                                                                              // Bjorken scaling slope
        Double_t E_beam = 6.5;                                                                              // TeV of energy for each proton
        Double_t mass = 2E-6;                                                                               // 2 MeV of mass assumed for the quark (in TeV)
        Double_t E_i1_CoM, E_i2_CoM;                                                                        // TeV of energy for the incoming partons (in CoM frame)
        Double_t sqrt_s;                                                                                    // √s 
        Double_t E1_prime_CoM, E2_prime_CoM;                                                                // TeV of energy for the outgoing partons (in CoM frame)
        Double_t E1_prime_LAB, E2_prime_LAB;                                                                // TeV of energy for the outgoing partons (in LAB frame)
        Double_t theta_CoM, phi1, phi2, theta1_LAB, theta2_LAB;                                             // angles of the outgoing partons
        TLorentzVector p1_prime, p2_prime, pSUM_CoM;                    // Lorentz vectors of partons after collision
        Double_t velocity_CoM_frame;                                                                        // velocity of CoM wrt to LAB
        
    
        // post-jets
        Double_t x1, x2, split_prob1, split_prob2;                                                          // Bjorken scalings and probability of jet splitting
        Double_t E_j1_CoM, E_j2_CoM, E_j3_CoM, E_j4_CoM;                                                    // TeV of energy of the jets (in CoM frame)
        Double_t E_j1_LAB, E_j2_LAB, E_j3_LAB, E_j4_LAB;                                                    // TeV of energy of the jets (in CoM frame)
        Double_t theta_j1_CoM, theta_j2_CoM, theta_j3_CoM, theta_j4_CoM, phi_j1_CoM, phi_j2_CoM, phi_j3_CoM, phi_j4_CoM;    // angles of the jets emitted in the CoM; j2 and j4 are collinear with j1 and j3 respectively
        Double_t theta_j1_LAB, theta_j2_LAB, theta_j3_LAB, theta_j4_LAB, phi_j1_LAB, phi_j2_LAB, phi_j3_LAB, phi_j4_LAB;    // angles of the jets emitted in the CoM; j2 and j4 are collinear with j1 and j3 respectively
        TLorentzVector p_j1, p_j2, p_j3, p_j4;                                              // Lorentz vectors of jets (in CoM frame)
        
        
        TF1 *F = new TF1("F", "0.5 + 0.5 * TMath::Erf((x-[0]) / [1])", 0., 1000.);                           // detection efficiency (domain in GeV)
        F->SetParameter(0, 55.);
        F->SetParameter(1, 10.);
        TF1 *B = new TF1("B", "0.7 * (0.5 + 2 * TMath::Erf((1.8 - TMath::Abs(x)) / 0.3))");                 // b-tagging efficiency
        Double_t F_j1, F_j2, F_j3, F_j4;                                                                    // uniformly random numbers to be compared to F distr (detection efficiency)
        Double_t F_th;                                                                                      // upper random value for all the F_j's
        Double_t B_j1, B_j2, B_j3, B_j4;                                                                    // uniformly random numbers to be compared to B distr (b-tagging efficiency)
        Double_t B_th;                                                                                      // upper random value for all the B_j's
        Int_t Nbtaggedjets;                                                                                 // nb of jets b-tagged


    // HISTOGRAMS
        TH1D *hist_x1 = new TH1D("", "Bjorken x1; x1; Events", 200, 0, 2);
        TH1D *hist_x2 = new TH1D("", "Bjorken x2; x2; Events", 200, 0, 2);
        

        TH1D *hist_E1_prime_CoM = new TH1D("", "E parton 1 (blue CoM and red LAB); TeV; Events", 250, -1, 20);
        hist_E1_prime_CoM->SetLineColor(kBlue);
        TH1D *hist_E1_prime_LAB = new TH1D("", "", 250, -1, 20);
        hist_E1_prime_LAB->SetLineColor(kRed);


        TH1D *hist_E2_prime_CoM = new TH1D("", "E parton 2 (blue CoM and red LAB); TeV; Events", 250, -1, 20);
        hist_E2_prime_CoM->SetLineColor(kBlue);
        TH1D *hist_E2_prime_LAB = new TH1D("", "", 250, -1, 20);
        hist_E2_prime_LAB->SetLineColor(kRed);


        TH1D *hist_sumE_CoM = new TH1D("", "sumE of 2partons (blue CoM and red LAB); TeV; Events", 250, -1, 40);
        hist_sumE_CoM->SetLineColor(kBlue);
        TH1D *hist_sumE_LAB = new TH1D("", "", 250, -1, 40);
        hist_sumE_LAB->SetLineColor(kRed);


        TH1D *hist_v = new TH1D("", "Velocity of the CoM frame (in c units); v/c; Events", 250, -1.2, 1.2);
        hist_v->SetLineColor(kBlue);


        TH1D *hist_thetaLAB1 = new TH1D("", "thetaLAB (blue parton1 and red parton2); #theta(rad); Events", 250, -2, TMath::Pi() + 1);
        hist_thetaLAB1->SetLineColor(kBlue);
        TH1D *hist_thetaLAB2 = new TH1D("", "", 250, -2, TMath::Pi() + 1);
        hist_thetaLAB2->SetLineColor(kRed);


        TH1D *hist_E4j_CoM = new TH1D("", "energy sum of 4jets (blue CoM and red LAB); TeV; Events", 250, -1., 40.);
        hist_E4j_CoM->SetLineColor(kBlue);
        TH1D *hist_E4j_LAB = new TH1D("", "", 250, -1., 40.);
        hist_E4j_LAB->SetLineColor(kRed);


        TH1D *hist_thetaj1_LAB = new TH1D("", "theta LAB 4jets; #theta(rad); Events", 250, -10., 10.);
        hist_thetaj1_LAB->SetLineColor(kBlue);
        TH1D *hist_thetaj2_LAB = new TH1D("", "", 250, -10., 10.);
        hist_thetaj2_LAB->SetLineColor(kRed);
        TH1D *hist_thetaj3_LAB = new TH1D("", "", 250, -10., 10.);
        hist_thetaj3_LAB->SetLineColor(kGreen);
        TH1D *hist_thetaj4_LAB = new TH1D("", "", 250, -10., 10.);
        hist_thetaj4_LAB->SetLineColor(kBlack);
        

        TH1D *hist_phij1_LAB = new TH1D("", "phi LAB 4jets; #phi(rad); Events", 250, -10., 10.);
        hist_phij1_LAB->SetLineColor(kBlue);
        TH1D *hist_phij2_LAB = new TH1D("", "", 250, -10., 10.);
        hist_phij2_LAB->SetLineColor(kRed);
        TH1D *hist_phij3_LAB = new TH1D("", "", 250, -10., 10.);
        hist_phij3_LAB->SetLineColor(kGreen);
        TH1D *hist_phij4_LAB = new TH1D("", "", 250, -10., 10.);
        hist_phij4_LAB->SetLineColor(kBlack);


        TH1D *hist_etaj1 = new TH1D("", "eta 4jets; #eta; Events", 250, -10., 10.);
        hist_etaj1->SetLineColor(kBlue);
        TH1D *hist_etaj2 = new TH1D("", "", 250, -10., 10.);
        hist_etaj2->SetLineColor(kRed);
        TH1D *hist_etaj3 = new TH1D("", "", 250, -10., 10.);
        hist_etaj3->SetLineColor(kGreen);
        TH1D *hist_etaj4 = new TH1D("", "", 250, -10., 10.);
        hist_etaj4->SetLineColor(kBlack);


        TH1D *hist_pTj1 = new TH1D("", "pT 4jets; p_{T} (TeV); Events", 250, 0, 10);
        hist_pTj1->SetLineColor(kBlue);
        TH1D *hist_pTj2 = new TH1D("", "", 250, 0, 10);
        hist_pTj2->SetLineColor(kRed);
        TH1D *hist_pTj3 = new TH1D("", "", 250, 0, 10);
        hist_pTj3->SetLineColor(kGreen);
        TH1D *hist_pTj4 = new TH1D("", "", 250, 0, 10);
        hist_pTj4->SetLineColor(kBlack);


        TH1D *hist_F = new TH1D("", "F", 250, 0., 1.);
        TH1D *hist_B = new TH1D("", "B", 250, 0., 1.);

    
    // create output file and output TTree
    TString outfname;
    outfname.Form("4jQCDbackground_autoBoost_%ievents.root", Nevents);
    TFile *outFile = new TFile(outfname, "recreate");
    TTree *Events_4jTree = new TTree("Events_4j", "Events_4j");
    
    jetEvent_t jetEvent1, jetEvent2, jetEvent3, jetEvent4;
    TBranch *j1 = Events_4jTree->Branch("j1", &jetEvent1);
    TBranch *j2 = Events_4jTree->Branch("j2", &jetEvent2);
    TBranch *j3 = Events_4jTree->Branch("j3", &jetEvent3);
    TBranch *j4 = Events_4jTree->Branch("j4", &jetEvent4);




    // event generation
    for (Int_t i = 0; i < Nevents; i++)
    {   
        x1 = gRandom->Exp(tau); // exp (-t / tau)
        x2 = gRandom->Exp(tau); // exp (-t / tau)
        while (x1 > 1. || x2 > 1.)
        {
            x1 = gRandom->Exp(tau);
            x2 = gRandom->Exp(tau);
        }

        // calculate energy of the incoming partons
        E_i1_CoM = x1 * E_beam;
        E_i2_CoM = x2 * E_beam;
        
        // energy in the CoM
        sqrt_s = E_i1_CoM + E_i2_CoM; 

        // equal energy partition between partons in the CoM
        E1_prime_CoM = sqrt_s / 2.; 
        E2_prime_CoM = sqrt_s / 2.;

        // angles of the outgoing partons in the 2->2 process
        theta_CoM = gRandom->Uniform(0, TMath::Pi());
        phi1 = gRandom->Uniform(0, 2 * TMath::Pi()); 
        phi2 = phi1 + TMath::Pi();

        // considering 0 mass also for the outgoing partons
        p1_prime.SetPxPyPzE(E1_prime_CoM * TMath::Sin(theta_CoM) * TMath::Cos(phi1), E1_prime_CoM * TMath::Sin(theta_CoM) * TMath::Sin(phi1), E1_prime_CoM * TMath::Cos(theta_CoM), E1_prime_CoM);
        p2_prime.SetPxPyPzE(E2_prime_CoM * TMath::Sin(TMath::Pi() - theta_CoM) * TMath::Cos(phi2), E2_prime_CoM * TMath::Sin(TMath::Pi() - theta_CoM) * TMath::Sin(phi2), E2_prime_CoM * TMath::Cos(TMath::Pi() - theta_CoM), E2_prime_CoM); // theta for second parton is (pi - theta) in CoM
        
        // 4-momentum of the composite system in the CoM - 3-momentum should be ~0
        pSUM_CoM.SetPxPyPzE(p1_prime.Px() + p2_prime.Px(), p2_prime.Py() + p2_prime.Py(), p1_prime.Pz() + p2_prime.Pz(), p1_prime.E() + p2_prime.E());
        // Indeed the CoM of the composite components are nearly 0 except rounding errors
        // cout << p1_prime_CoM.Px() + p2_prime_CoM.Px() << " " << p1_prime_CoM.Py() + p2_prime_CoM.Py() << " " << p1_prime_CoM.Pz() + p2_prime_CoM.Pz() << " " << p1_prime_CoM.E() + p2_prime_CoM.E() << endl;
        

        // boosting momenta to the LAB frame
        p1_prime.Boost(0., 0., x1 - x2);
        p2_prime.Boost(0., 0., x1 - x2);

        // calculating energies in the LAB frame (just for visualization, they are implicitly calculated in the boost_toLAB)
        E1_prime_LAB = p1_prime.E();
        E2_prime_LAB = p2_prime.E();


        theta1_LAB = p1_prime.Theta();
        theta2_LAB = p2_prime.Theta();





        // histograms up to before jets generation
        hist_x1->Fill(x1);
        hist_x2->Fill(x2);

        hist_E1_prime_CoM->Fill(E1_prime_CoM);
        hist_E1_prime_LAB->Fill(E1_prime_LAB);
        hist_E2_prime_CoM->Fill(E2_prime_CoM);
        hist_E2_prime_LAB->Fill(E2_prime_LAB);

        hist_sumE_CoM->Fill(E1_prime_CoM + E2_prime_CoM);
        hist_sumE_LAB->Fill(E1_prime_LAB + E2_prime_LAB);

        hist_v->Fill(x1 - x2);

        hist_thetaLAB1->Fill(theta1_LAB);
        hist_thetaLAB2->Fill(theta2_LAB);






        // probabilities of splitting
        split_prob1 = gRandom->Uniform(-1, 1);
        split_prob2 = gRandom->Uniform(-1, 1);

        if (split_prob1 > 0 && split_prob2 > 0) // splitting to 4 partons
        {

            // will use split prob also as a separator of energy since in the case of splitting split_prob is in (0, 1)
            E_j1_CoM = split_prob1 * E1_prime_CoM; // using the energy in CoM !!
            E_j2_CoM = (1 - split_prob1) * E1_prime_CoM;
            E_j3_CoM = split_prob2 * E2_prime_CoM;
            E_j4_CoM = (1 - split_prob2) * E2_prime_CoM;

            hist_E4j_CoM->Fill(E_j1_CoM + E_j2_CoM + E_j3_CoM + E_j4_CoM);

            // theta angles of jets in the CoM
            theta_j1_CoM = gRandom->Uniform(0, TMath::Pi());
            theta_j2_CoM = TMath::Pi() - theta_j1_CoM;
            theta_j3_CoM = gRandom->Uniform(0, TMath::Pi());
            theta_j4_CoM = TMath::Pi() - theta_j3_CoM;

            // phi angles of jets in the CoM
            phi_j1_CoM = gRandom->Uniform(0, 2 * TMath::Pi());
            phi_j2_CoM = phi_j1_CoM + TMath::Pi();
            phi_j3_CoM = gRandom->Uniform(0, 2 * TMath::Pi());
            phi_j4_CoM = phi_j3_CoM + TMath::Pi();

            // CoM 4-momenta of jets
            p_j1.SetPxPyPzE(E_j1_CoM * TMath::Sin(theta_j1_CoM) * TMath::Cos(phi_j1_CoM), E_j1_CoM * TMath::Sin(theta_j1_CoM) * TMath::Sin(phi_j1_CoM), E_j1_CoM * TMath::Cos(theta_j1_CoM), E_j1_CoM);
            p_j2.SetPxPyPzE(E_j2_CoM * TMath::Sin(theta_j2_CoM) * TMath::Cos(phi_j2_CoM), E_j2_CoM * TMath::Sin(theta_j2_CoM) * TMath::Sin(phi_j2_CoM), E_j2_CoM * TMath::Cos(theta_j2_CoM), E_j2_CoM);
            p_j3.SetPxPyPzE(E_j3_CoM * TMath::Sin(theta_j3_CoM) * TMath::Cos(phi_j3_CoM), E_j3_CoM * TMath::Sin(theta_j3_CoM) * TMath::Sin(phi_j3_CoM), E_j3_CoM * TMath::Cos(theta_j3_CoM), E_j3_CoM);
            p_j4.SetPxPyPzE(E_j4_CoM * TMath::Sin(theta_j4_CoM) * TMath::Cos(phi_j4_CoM), E_j4_CoM * TMath::Sin(theta_j4_CoM) * TMath::Sin(phi_j4_CoM), E_j4_CoM * TMath::Cos(theta_j4_CoM), E_j4_CoM);

            
            // boosting jets 4-momenta to the LAB frame
            p_j1.Boost(0., 0., x1 - x2);
            p_j2.Boost(0., 0., x1 - x2);
            p_j3.Boost(0., 0., x1 - x2);
            p_j4.Boost(0., 0., x1 - x2);

            
            // thetaLAB of jets
            theta_j1_LAB = p_j1.Theta();
            theta_j2_LAB = p_j2.Theta();
            theta_j3_LAB = p_j3.Theta();
            theta_j4_LAB = p_j4.Theta();

            hist_thetaj1_LAB->Fill(theta_j1_LAB);
            hist_thetaj2_LAB->Fill(theta_j2_LAB);
            hist_thetaj3_LAB->Fill(theta_j3_LAB);
            hist_thetaj4_LAB->Fill(theta_j4_LAB);


            // pT and eta
            Double_t pT_j1, pT_j2, pT_j3, pT_j4;
            Double_t eta_j1, eta_j2, eta_j3, eta_j4;

            pT_j1 = TMath::Sqrt(p_j1.Px()*p_j1.Px() + p_j1.Py()*p_j1.Py());
            pT_j2 = TMath::Sqrt(p_j2.Px()*p_j2.Px() + p_j2.Py()*p_j2.Py());
            pT_j3 = TMath::Sqrt(p_j3.Px()*p_j3.Px() + p_j3.Py()*p_j3.Py());
            pT_j4 = TMath::Sqrt(p_j4.Px()*p_j4.Px() + p_j4.Py()*p_j4.Py());

            hist_pTj1->Fill(pT_j1);
            hist_pTj2->Fill(pT_j2);
            hist_pTj3->Fill(pT_j3);
            hist_pTj4->Fill(pT_j4);

            eta_j1 = - TMath::Log(TMath::Tan(theta_j1_LAB / 2));
            eta_j2 = - TMath::Log(TMath::Tan(theta_j2_LAB / 2));
            eta_j3 = - TMath::Log(TMath::Tan(theta_j3_LAB / 2));
            eta_j4 = - TMath::Log(TMath::Tan(theta_j4_LAB / 2));

            hist_etaj1->Fill(eta_j1);
            hist_etaj2->Fill(eta_j2);
            hist_etaj3->Fill(eta_j3);
            hist_etaj4->Fill(eta_j4);

            phi_j1_LAB = 2 * TMath::ATan(p_j1.Py() / p_j1.Px());
            phi_j2_LAB = 2 * TMath::ATan(p_j2.Py() / p_j2.Px());
            phi_j3_LAB = 2 * TMath::ATan(p_j3.Py() / p_j3.Px());
            phi_j4_LAB = 2 * TMath::ATan(p_j4.Py() / p_j4.Px());

            hist_phij1_LAB->Fill(phi_j1_LAB);
            hist_phij2_LAB->Fill(phi_j2_LAB);
            hist_phij3_LAB->Fill(phi_j3_LAB);
            hist_phij4_LAB->Fill(phi_j4_LAB);

            E_j1_LAB = p_j1.E();
            E_j2_LAB = p_j2.E();
            E_j3_LAB = p_j3.E();
            E_j4_LAB = p_j4.E();

            hist_E4j_LAB->Fill(E_j1_LAB + E_j2_LAB + E_j3_LAB + E_j4_LAB);

            p_j1.SetPtEtaPhiE(pT_j1, eta_j1, phi_j1_LAB, E_j1_LAB);
            p_j2.SetPtEtaPhiE(pT_j2, eta_j2, phi_j2_LAB, E_j2_LAB);
            p_j3.SetPtEtaPhiE(pT_j3, eta_j3, phi_j3_LAB, E_j3_LAB);
            p_j4.SetPtEtaPhiE(pT_j4, eta_j4, phi_j4_LAB, E_j4_LAB);



            // detection acceptance
            F_j1 = gRandom->Uniform(0, 1);
            F_j2 = gRandom->Uniform(0, 1);
            F_j3 = gRandom->Uniform(0, 1);
            F_j4 = gRandom->Uniform(0, 1);

            if (F_j1 < F->Eval(p_j1.Pt()*1E3) && F_j2 < F->Eval(p_j2.Pt()*1E3) && F_j3 < F->Eval(p_j3.Pt()*1E3) && F_j4 < F->Eval(p_j4.Pt()*1E3)) // all 4 jets detected
            {
                // b-tagging efficiency
                B_j1 = gRandom->Uniform(0, 1);
                B_j2 = gRandom->Uniform(0, 1);
                B_j3 = gRandom->Uniform(0, 1);
                B_j4 = gRandom->Uniform(0, 1);
                B_th = F->GetRandom();

                if (B_j1 < B_th && B_j2 < B_th && B_j3 < B_th && B_j4 < B_th) // all 4 jets b-tagged
                {   
                    Nbtaggedjets = 4;
                    jetEvent1.pt = p_j1.Pt();
                    jetEvent1.eta = p_j1.Eta();
                    jetEvent1.phi = p_j1.Phi();

                    jetEvent2.pt = p_j2.Pt();
                    jetEvent2.eta = p_j2.Eta();
                    jetEvent2.phi = p_j2.Phi();

                    jetEvent3.pt = p_j3.Pt();
                    jetEvent3.eta = p_j3.Eta();
                    jetEvent3.phi = p_j3.Phi();

                    jetEvent4.pt = p_j4.Pt();
                    jetEvent4.eta = p_j4.Eta();
                    jetEvent4.phi = p_j4.Phi();
                    
                    
                    Events_4jTree->Fill();
                }

                if ((B_j1 < B_th && B_j2 < B_th && B_j3 < B_th && B_j4 > B_th) || (B_j1 < B_th && B_j2 < B_th && B_j3 > B_th && B_j4 < B_th) || (B_j1 < B_th && B_j2 > B_th && B_j3 < B_th && B_j4 < B_th) || ((B_j1 > B_th && B_j2 < B_th && B_j3 < B_th && B_j4 < B_th))) // 3 jets b-tagged
                {
                    Nbtaggedjets = 3;
                    //mytuple->Fill(p_j1_CoM.Px(), p_j1_CoM.Py(), p_j1_CoM.Pz(), p_j1_CoM.E(), p_j1_CoM.Eta(), p_j1_CoM.Phi(), p_j2_CoM.Px(), p_j2_CoM.Py(), p_j2_CoM.Pz(), p_j2_CoM.E(), p_j2_CoM.Eta(), p_j2_CoM.Phi(), p_j3_CoM.Px(), p_j3_CoM.Py(), p_j3_CoM.Pz(), p_j3_CoM.E(), p_j3_CoM.Eta(), p_j3_CoM.Phi(), p_j4_CoM.Px(), p_j4_CoM.Py(), p_j4_CoM.Pz(), p_j4_CoM.E(), p_j4_CoM.Eta(), p_j4_CoM.Phi(), Nbtaggedjets)
                }

            }
            else
            {
                continue;
            }

        }
        else if ((split_prob1 > 0 && split_prob2 < 0) || (split_prob1 < 0 && split_prob2 > 0)) // splitting to three partons
        {
            //cout << "three jets" << endl;
        }
        else if (split_prob1 < 0 && split_prob2 < 0) // no splitting (two partons)
        {
            //cout << "two jets" << endl;
        }
    }
    outFile->cd();
    Events_4jTree->Write("", TObject::kOverwrite);
    outFile->Close();
    delete outFile;

    // plotting for checking
    TCanvas *c1 = new TCanvas("c1", "x1");
    c1->SetLogy();
    hist_x1->Draw();

    TCanvas *c2 = new TCanvas("c2", "x2");
    c2->SetLogy();
    hist_x2->Draw();


    TCanvas *c6 = new TCanvas("c6", "Vel. CoM");
    c6->SetLogy();
    hist_v->Draw();


    TCanvas *c3 = new TCanvas("c3", "E1");
    c3->SetLogy();
    hist_E1_prime_CoM->Draw();
    hist_E1_prime_LAB->Draw("same");

    TCanvas *c4 = new TCanvas("c4", "E2");
    c4->SetLogy();
    hist_E2_prime_CoM->Draw();
    hist_E2_prime_LAB->Draw("same");


    TCanvas *c7 = new TCanvas("c7", "theta_LAB");
    hist_thetaLAB1->Draw();
    hist_thetaLAB2->Draw("same");


    TCanvas *c8 = new TCanvas("c8", "sumE 2partons");
    c8->SetLogy();
    hist_sumE_CoM->Draw();
    hist_sumE_LAB->Draw("same");


    TCanvas *c9 = new TCanvas("c9", "sumE 4jets");
    c9->SetLogy();
    hist_E4j_CoM->Draw("");
    hist_E4j_LAB->Draw("same");


    TCanvas *c10 = new TCanvas("c10", "thetaLAB 4jets");
    hist_thetaj1_LAB->Draw("");
    hist_thetaj2_LAB->Draw("same");
    hist_thetaj3_LAB->Draw("same");
    hist_thetaj4_LAB->Draw("same");

    TCanvas *c11 = new TCanvas("c11", "phiLAB 4jets");
    hist_phij1_LAB->Draw("");
    hist_phij2_LAB->Draw("same");
    hist_phij3_LAB->Draw("same");
    hist_phij4_LAB->Draw("same");

    TCanvas *c12 = new TCanvas("c12", "eta 4jets");
    hist_etaj1->Draw("");
    hist_etaj2->Draw("same");
    hist_etaj3->Draw("same");
    hist_etaj4->Draw("same");

    TCanvas *c14 = new TCanvas("c14", "F");
    hist_F->Draw();

    TCanvas *c13 = new TCanvas("c13", "pT 4jets");
    c12->SetLogy();
    hist_pTj1->Draw("");
    hist_pTj2->Draw("same");
    hist_pTj3->Draw("same");
    hist_pTj4->Draw("same");
}