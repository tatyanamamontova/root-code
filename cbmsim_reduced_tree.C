#define cbmsim_reduced_tree_cxx
#include "cbmsim_reduced_tree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>

void cbmsim_reduced_tree::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L cbmsim_reduced_tree.C
//      Root > cbmsim_reduced_tree t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) 
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      //! loop for mpd_reconstruction
      for (Long64_t itracks=0;itracks<n_tracks_mpd;itracks++)
      {
        //! common cuts
      	if(TMath::Abs(signed_pt_mpd[itracks])<=3  && TMath::Abs(eta_mpd[itracks]) <= 1.5  && n_hits_mpd[itracks] >= 32 && mother_ID_mc[id_from_mc_mpd[itracks]] == -1) 
      	{
          //! for protons
          if (PDG_code_mc[id_from_mc_mpd[itracks]]==2212)
          {
            h_pt_mpd_p->Fill(TMath::Abs(signed_pt_mpd[itracks]));
	    h_eta_mpd_p->Fill(eta_mpd[itracks]);
	    h_corr_pt_p->Fill(TMath::Abs(signed_pt_mpd[itracks]),
		(TMath::Abs(signed_pt_mpd[itracks]) - pt_mc[id_from_mc_mpd[itracks]])/(TMath::Abs(signed_pt_mpd[itracks]) + pt_mc[id_from_mc_mpd[itracks]]));
	    h_prof_p->Fill(TMath::Abs(signed_pt_mpd[itracks]),
		(TMath::Abs(signed_pt_mpd[itracks]) - pt_mc[id_from_mc_mpd[itracks]])/(TMath::Abs(signed_pt_mpd[itracks]) + pt_mc[id_from_mc_mpd[itracks]]),1);

          }
          //! for kaons
          if (PDG_code_mc[id_from_mc_mpd[itracks]]==321)
          {
            h_pt_mpd_k->Fill(TMath::Abs(signed_pt_mpd[itracks]));
	    h_eta_mpd_k->Fill(eta_mpd[itracks]);
	    h_corr_pt_k->Fill(TMath::Abs(signed_pt_mpd[itracks]),
		(TMath::Abs(signed_pt_mpd[itracks]) - pt_mc[id_from_mc_mpd[itracks]])/(TMath::Abs(signed_pt_mpd[itracks]) + pt_mc[id_from_mc_mpd[itracks]]));
            h_prof_k->Fill(TMath::Abs(signed_pt_mpd[itracks]),
		(TMath::Abs(signed_pt_mpd[itracks]) - pt_mc[id_from_mc_mpd[itracks]])/(TMath::Abs(signed_pt_mpd[itracks]) + pt_mc[id_from_mc_mpd[itracks]]),1);
          }
          //! for pions
          if (PDG_code_mc[id_from_mc_mpd[itracks]]==211)
          {
            h_pt_mpd_pi->Fill(TMath::Abs(signed_pt_mpd[itracks]));
	    h_eta_mpd_pi->Fill(eta_mpd[itracks]);
            h_corr_pt_pi->Fill(TMath::Abs(signed_pt_mpd[itracks]),
		(TMath::Abs(signed_pt_mpd[itracks]) - pt_mc[id_from_mc_mpd[itracks]])/(TMath::Abs(signed_pt_mpd[itracks]) + pt_mc[id_from_mc_mpd[itracks]]));
            h_prof_pi->Fill(TMath::Abs(signed_pt_mpd[itracks]),
		(TMath::Abs(signed_pt_mpd[itracks]) - pt_mc[id_from_mc_mpd[itracks]])/(TMath::Abs(signed_pt_mpd[itracks]) + pt_mc[id_from_mc_mpd[itracks]]),1);
          }
      		
   	}

      }
      //! the same for mc-simulations  
      for (Long64_t itracks=0;itracks<n_tracks_mc;itracks++)
      {
        if(pt_mc[itracks]<=3 && TMath::Abs(eta_mc[itracks]) <= 1.5  && mother_ID_mc[itracks] == -1) 
       	{
          if(PDG_code_mc[itracks]==2212)
          {
            h_pt_mc_p->Fill(pt_mc[itracks]);
	    h_eta_mc_p->Fill(eta_mc[itracks]);
          } 

          if(PDG_code_mc[itracks]==321)
          {
            h_pt_mc_k->Fill(pt_mc[itracks]);
	    h_eta_mc_k->Fill(eta_mc[itracks]);
          }

          if(PDG_code_mc[itracks]==211)
          {
            h_pt_mc_pi->Fill(pt_mc[itracks]);
	    h_eta_mc_pi->Fill(eta_mc[itracks]);
          } 

        }

      }

    }
   h_eff_p->Divide(h_pt_mpd_p, h_pt_mc_p);
   h_eff_k->Divide(h_pt_mpd_k, h_pt_mc_k);
   h_eff_pi->Divide(h_pt_mpd_pi, h_pt_mc_pi);

   h_eta_eff_p->Divide(h_eta_mpd_p, h_eta_mc_p);
   h_eta_eff_k->Divide(h_eta_mpd_k, h_eta_mc_k);
   h_eta_eff_pi->Divide(h_eta_mpd_pi, h_eta_mc_pi);


   //! for tranverse momentum
   h_eff_p->SetMarkerStyle(kDot);
   h_eff_k->SetMarkerStyle(kCircle);
   h_eff_pi->SetMarkerStyle(kStar);
   h_eff_p->SetMarkerColor(kRed);
   h_eff_k->SetMarkerColor(kGreen);
   h_eff_pi->SetMarkerColor(kBlue);
   h_eff_p->SetLineColor(kRed);
   h_eff_k->SetLineColor(kGreen);
   h_eff_pi->SetLineColor(kBlue);

  	
   TCanvas *c1 = new TCanvas("efficiency of p_T","efficiency of p_T", 600, 500);
   TH1F *frame1 = (TH1F*)c1->DrawFrame(0,-1,3,7);
   frame1->GetXaxis()->SetTitle("p_T [GeV]");
   frame1->GetYaxis()->SetTitle("efficiency");
   h_eff_p->Draw("same");
   h_eff_k->Draw("same");
   h_eff_pi->Draw("same");

   TLegend *legend1 = new TLegend(0.55,0.65,0.16,0.82);
   legend1->SetHeader("Efficiencies of particle reconstruction");
   legend1->AddEntry(h_eff_p,"protons","lep");
   legend1->AddEntry(h_eff_k,"kaons","lep");
   legend1->AddEntry(h_eff_pi,"pions","lep");
   legend1->Draw();

   //! For pseudorapidity

   h_eta_eff_p->SetMarkerStyle(kDot);
   h_eta_eff_k->SetMarkerStyle(kCircle);
   h_eta_eff_pi->SetMarkerStyle(kStar);
   h_eta_eff_p->SetMarkerColor(kRed);
   h_eta_eff_k->SetMarkerColor(kGreen);
   h_eta_eff_pi->SetMarkerColor(kBlue);
   h_eta_eff_p->SetLineColor(kRed);
   h_eta_eff_k->SetLineColor(kGreen);
   h_eta_eff_pi->SetLineColor(kBlue);

  	
   TCanvas *c2 = new TCanvas("efficiency of pseudorapidity","efficiency of pseudorapidity", 600, 500);
   TH1F *frame2 = (TH1F*)c2->DrawFrame(0,-1,3,7);
   frame2->GetXaxis()->SetTitle("eta");
   frame2->GetYaxis()->SetTitle("efficiency");
   h_eta_eff_p->Draw("same");
   h_eta_eff_k->Draw("same");
   h_eta_eff_pi->Draw("same");

   TLegend *legend2 = new TLegend(0.55,0.65,0.16,0.82);
   legend2->SetHeader("Efficiencies of particle reconstruction");
   legend2->AddEntry(h_eff_p,"protons","lep");
   legend2->AddEntry(h_eff_k,"kaons","lep");
   legend2->AddEntry(h_eff_pi,"pions","lep");
   legend2->Draw();
   

//! For correlation of tranverse momentum
   TCanvas *c3 = new TCanvas("correlation of tranverse momentum for proton","correlation of tranverse momentum for proton", 600, 500);	
   h_corr_pt_p->Draw("colz");
   TCanvas *c4 = new TCanvas("correlation of tranverse momentum for kaon","correlation of tranverse momentum for kaon", 600, 500);
   h_corr_pt_k->Draw("colz");
   TCanvas *c5 = new TCanvas("correlation of tranverse momentum for pion","correlation of tranverse momentum for pion", 600, 500);
   h_corr_pt_pi->Draw("colz");


//! For profile histos
   TCanvas *c6 = new TCanvas("profile for correlation for proton","profile for correlation for proton", 200, 10, 600, 500);	
   h_prof_p->Draw();
   TCanvas *c7 = new TCanvas("profile for correlation for kaon","profile for correlation for kaon",  200, 10, 600, 500);
   h_prof_k->Draw();
   TCanvas *c8 = new TCanvas("profile for correlation for pion","profile for correlation for pion",  200, 10, 600, 500);
   h_prof_pi->Draw();

   TFile* output = new TFile("output.root", "recreate");
   output->cd();
   h_pt_mc_p->Write();
   h_pt_mpd_p->Write();
   h_eff_p->Write();
   h_pt_mc_k->Write();
   h_pt_mpd_k->Write();
   h_eff_k->Write();
   h_pt_mc_pi->Write();
   h_pt_mpd_pi->Write();
   h_eff_pi->Write();
   h_eta_mc_p->Write();
   h_eta_mpd_p->Write();
   h_eta_eff_p->Write(); 
   h_eta_mc_k->Write();
   h_eta_mpd_k->Write(); 
   h_eta_eff_k->Write(); 
   h_eta_mc_pi->Write(); 
   h_eta_mpd_pi->Write(); 
   h_eta_eff_pi->Write();
   h_corr_pt_p->Write(); 
   h_corr_pt_k->Write(); 
   h_corr_pt_pi->Write();
   h_prof_p->Write();
   h_prof_k->Write();
   h_prof_pi->Write(); 
   output->Close();

}
