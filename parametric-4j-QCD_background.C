/////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Javier Mariño Villadamigo                                                           //
//                                                                                             //
// Date: 30/11/2022                                                                            //
//                                                                                             //
// Thesis: parametric model for the generation of a 4-jets QCD background                      //
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
// NOTE:
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////



void parametric_4j_QCD_background() {

    delete gRandom;
    //int seed = 0;
    gRandom = new TRandom2();

    

    TNtuple *mytuple = new TNtuple("mytuple", "parametric_4j_QCD_background", "j_Px:j_Py:j_Pz:j_E:j_eta:j_phi");
    TFile *outputFile = new TFile("NTuple.root", "recreate");

    // VARIABLES

    Double_t tau = 1./10.; // Bjorken scaling slope
    Double_t E1 = 6.5; // TeV of energy for each proton
    Double_t E1_prime, E2_prime; // TeV of energy for the outgoing partons
    Double_t E_j1, E_j2, E_j3, E_j4; // TeV of energy of the jets

    Double_t x1, x2, y, split_prob1, split_prob2;
    Double_t q0;
    Double_t theta_CoM, phi1, phi2, theta1_LAB, theta2_LAB; // angles of the outgoing partons
    Double_t theta_j1, theta_j2, theta_j3, theta_j4,  phi_j1, phi_j2, phi_j3, phi_j4; // angles of the jets emitted in the CoM; j2 and j4 are collinear with j1 and j3 respectively
    
    Double_t F_j1, F_j2, F_j3, F_j4; // uniformly random numbers to be compared to F distr (detection efficiency)
    Double_t F_th; // upper random value for all the F_j's

    Double_t B_j1, B_j2, B_j3, B_j4; // uniformly random numbers to be compared to B distr (b-tagging efficiency)
    Double_t B_th; // upper random value for all the B_j's

    // detection efficiency
    TF1 *F = new TF1("F", "0.5 + 2 * (TMath::Erf(50-x) + TMath::Erf(30))");

    // b-tagging efficiency
    TF1 *B = new TF1("B", "0.7 * (0.5 + 2 * TMath::Erf((1.8 - TMath::Abs(x)) / 0.3))");

    Int_t Nbtaggedjets;

    // Lorentz vectors of partons before and after collision
    TLorentzVector p1(0, 0, E1, E1);
    TLorentzVector p2(0, 0, -E1, E1);
    TLorentzVector p1_prime, p2_prime;

    // Lorentz vectors of jets
    TLorentzVector p_j1, p_j2, p_j3, p_j4;


    TH1D *hist_x1 = new TH1D("", "Bjorken x1", 200, 0, 1);
    TH1D *hist_x2 = new TH1D("", "Bjorken x2", 200, 0, 1);
    TH1D *hist_E1_prime = new TH1D("", "E' parton 1", 250, 0, 20);
    TH1D *hist_E2_prime = new TH1D("", "E' parton 2", 250, 0, 20);
    TH1D *hist_sum = new TH1D("", "sum(E')", 250, 0, 20);
    TH1D *hist_boost = new TH1D("", "boost (along z) x1-x2", 250, -1.2, 1.2);

    TH1D *hist_thetaLAB1 = new TH1D("", "thetaLAB for both partons", 250, -1, TMath::Pi() + 1);
    hist_thetaLAB1->SetLineColor(kBlue);
    TH1D *hist_thetaLAB2 = new TH1D("", "thetaLAB", 250, -1, TMath::Pi() + 1);
    hist_thetaLAB2->SetLineColor(kRed);

    TH1D *hist_E4j = new TH1D("", "energy sum of 4j", 250, 0., 20.);

    TH1D *hist_F = new TH1D("", "F", 250, 0., 1.);
    TH1D *hist_B = new TH1D("", "B", 250, 0., 1.);



    for (Int_t i = 0; i<100000; i++)
    {
        x1 = gRandom->Exp(tau); // exp (t / tau)
        x2 = gRandom->Exp(tau); // exp (t / tau)
        y = gRandom->Uniform(-1, 1); // fraction of momentum transferred by one parton to another (in [-1, 1] to allow bidirectionality)

        q0 = y * E1; // energy transferred from one parton to another

        E1_prime = E1 - q0;
        E2_prime = E1 + q0;
        theta_CoM = gRandom->Uniform(0, TMath::Pi());
        phi1 = gRandom->Uniform(0, TMath::Pi()); // azimuthal angle for one of the two partons in the 2->2 hard process
        phi2 = phi1 + TMath::Pi();

        // considering 0 mass also for the outgoing partons
        p1_prime.SetPxPyPzE(E1_prime * TMath::Sin(theta_CoM) * TMath::Cos(phi1), E1_prime * TMath::Sin(theta_CoM) * TMath::Sin(phi1), E1_prime * TMath::Cos(theta_CoM), E1_prime);
        p2_prime.SetPxPyPzE(E2_prime * TMath::Sin(TMath::Pi() - theta_CoM) * TMath::Cos(phi2), E2_prime * TMath::Sin(TMath::Pi() - theta_CoM) * TMath::Sin(phi2), E2_prime * TMath::Cos(TMath::Pi() - theta_CoM), E2_prime); // theta for second parton is (pi - theta) in CoM




        /////
        p1_prime.Boost(0, 0, x1-x2); // *********** SOMETHING WRONG WITH THE BOOST *********+
        p2_prime.Boost(0, 0, x1-x2);
        /////

        theta1_LAB = p1_prime.Theta();
        theta2_LAB = p2_prime.Theta();

        split_prob1 = gRandom->Uniform(-1, 1);
        split_prob2 = gRandom->Uniform(-1, 1);

        // histograms up to before jets generation
        hist_x1->Fill(x1);
        hist_x2->Fill(x2);
        hist_E1_prime->Fill(E1_prime);
        hist_E2_prime->Fill(E2_prime);
        hist_sum->Fill(E1_prime + E2_prime);
        hist_boost->Fill(x1-x2);
        hist_thetaLAB1->Fill(theta1_LAB);
        hist_thetaLAB2->Fill(theta2_LAB);


        // probabilities of splitting
        if (split_prob1 > 0 && split_prob2 > 0) // splitting to 4 partons
        {
            //cout << "four jets" << endl;

            // will use split prob also as a separator of energy
            E_j1 = split_prob1 * E1_prime; // using the energy in CoM !!
            E_j2 = (1 - split_prob1) * E1_prime;
            E_j3 = split_prob2 * E2_prime;
            E_j4 = (1 - split_prob2) * E2_prime;

            // theta angles of jets in the CoM
            theta_j1 = gRandom->Uniform(0, TMath::Pi());
            theta_j2 = TMath::Pi() - theta_j1;
            theta_j3 = gRandom->Uniform(0, TMath::Pi());
            theta_j4 = TMath::Pi() - theta_j3;

            // phi angles of jets in the CoM
            phi_j1 = gRandom->Uniform(0, TMath::Pi());
            phi_j2 = phi_j1 + TMath::Pi();
            phi_j3 = gRandom->Uniform(0, TMath::Pi());
            phi_j4 = phi_j3 + TMath::Pi();

            // 4-momenta of jets
            p_j1.SetPxPyPzE(E_j1 * TMath::Sin(theta_j1) * TMath::Cos(phi_j1), E_j1 * TMath::Sin(theta_j1) * TMath::Sin(phi_j1), E_j1 * TMath::Cos(theta_j1), E_j1);
            p_j2.SetPxPyPzE(E_j2 * TMath::Sin(theta_j2) * TMath::Cos(phi_j2), E_j2 * TMath::Sin(theta_j2) * TMath::Sin(phi_j2), E_j2 * TMath::Cos(theta_j2), E_j2);
            p_j3.SetPxPyPzE(E_j3 * TMath::Sin(theta_j3) * TMath::Cos(phi_j3), E_j3 * TMath::Sin(theta_j3) * TMath::Sin(phi_j3), E_j3 * TMath::Cos(theta_j3), E_j3);
            p_j4.SetPxPyPzE(E_j4 * TMath::Sin(theta_j4) * TMath::Cos(phi_j1), E_j4 * TMath::Sin(theta_j4) * TMath::Sin(phi_j4), E_j4 * TMath::Cos(theta_j4), E_j4);

            
            // checking conservation of energy in the CoM
            hist_E4j->Fill(p_j1.E() + p_j2.E() + p_j3.E() + p_j4.E());

            // boosting the jets
            p_j1.Boost(0, 0, x1-x2);
            p_j2.Boost(0, 0, x1-x2);
            p_j3.Boost(0, 0, x1-x2);
            p_j4.Boost(0, 0, x1-x2);

            F_j1 = gRandom->Uniform(0, 1);
            F_j2 = gRandom->Uniform(0, 1);
            F_j3 = gRandom->Uniform(0, 1);
            F_j4 = gRandom->Uniform(0, 1);
            F_th = F->GetRandom();
            hist_F->Fill(F_th);

            if (F_j1 < F_th && F_j2 < F_th && F_j3 < F_th && F_j4 < F_th) // all 4 jets detected
            {
                B_j1 = gRandom->Uniform(0, 1);
                B_j2 = gRandom->Uniform(0, 1);
                B_j3 = gRandom->Uniform(0, 1);
                B_j4 = gRandom->Uniform(0, 1);
                B_th = F->GetRandom();

                if (B_j1 < B_th && B_j2 < B_th && B_j3 < B_th && B_j4 < B_th) // all 4 jets b-tagged
                {   
                    Nbtaggedjets = 4;
                    mytuple->Fill(p_j1.Px(), p_j1.Py(), p_j1.Pz(), p_j1.E(), p_j1.Eta(), p_j1.Phi());
                    mytuple->Fill(p_j2.Px(), p_j2.Py(), p_j2.Pz(), p_j2.E(), p_j2.Eta(), p_j2.Phi());
                    mytuple->Fill(p_j3.Px(), p_j3.Py(), p_j3.Pz(), p_j3.E(), p_j3.Eta(), p_j3.Phi());
                    mytuple->Fill(p_j4.Px(), p_j4.Py(), p_j4.Pz(), p_j4.E(), p_j4.Eta(), p_j4.Phi());
                }

                if ((B_j1 < B_th && B_j2 < B_th && B_j3 < B_th && B_j4 > B_th) || (B_j1 < B_th && B_j2 < B_th && B_j3 > B_th && B_j4 < B_th) || (B_j1 < B_th && B_j2 > B_th && B_j3 < B_th && B_j4 < B_th) || ((B_j1 > B_th && B_j2 < B_th && B_j3 < B_th && B_j4 < B_th))) // 3 jets b-tagged
                {
                    Nbtaggedjets = 3;
                    //mytuple->Fill(p_j1.Px(), p_j1.Py(), p_j1.Pz(), p_j1.E(), p_j1.Eta(), p_j1.Phi(), p_j2.Px(), p_j2.Py(), p_j2.Pz(), p_j2.E(), p_j2.Eta(), p_j2.Phi(), p_j3.Px(), p_j3.Py(), p_j3.Pz(), p_j3.E(), p_j3.Eta(), p_j3.Phi(), p_j4.Px(), p_j4.Py(), p_j4.Pz(), p_j4.E(), p_j4.Eta(), p_j4.Phi(), Nbtaggedjets)
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











    // plotting for checking
    TCanvas *c1 = new TCanvas("c1", "x1");
    hist_x1->Draw();
    TCanvas *c2 = new TCanvas("c2", "x2");
    hist_x2->Draw();
    TCanvas *c3 = new TCanvas("c3", "E1_prime");
    hist_E1_prime->Draw();
    TCanvas *c4 = new TCanvas("c4", "E2_prime");
    hist_E2_prime->Draw();
    TCanvas *c5 = new TCanvas("c5", "E1_prime+E2_prime");
    hist_sum->Draw();
    TCanvas *c6 = new TCanvas("c6", "boost");
    c6->SetLogy();
    hist_boost->GetXaxis()->SetTitle("x1 - x2");
    hist_boost->Draw();
    TCanvas *c7 = new TCanvas("c7", "theta_LAB");
    hist_thetaLAB1->Draw();
    hist_thetaLAB2->Draw("same");
    TCanvas *c8 = new TCanvas("c8", "E4j");
    hist_E4j->Draw();
    TCanvas *c9 = new TCanvas("c9", "F");
    hist_F->Draw();
    TCanvas *c10 = new TCanvas("c10", "B");
    hist_B->Draw();
    
    
    mytuple->Write();
    outputFile->Close();
    













    //Double_t Q2 = x1 * x2 * s;
}