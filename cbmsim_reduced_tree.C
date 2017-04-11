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
      	if(TMath::Abs(signed_pt_mpd[itracks])<=3  && TMath::Abs(eta_mpd[itracks]) < 1.5  && n_hits_mpd[itracks] > 32 && mother_ID_mc[id_from_mc_mpd[itracks]] == -1) 
      	{
          //! for protons
          if (PDG_code_mc[id_from_mc_mpd[itracks]]==2212)
          {
            h_pt_mpd_p->Fill(TMath::Abs(signed_pt_mpd[itracks]));
          }
          //! for kaons
          if (PDG_code_mc[id_from_mc_mpd[itracks]]==130)
          {
            h_pt_mpd_k->Fill(TMath::Abs(signed_pt_mpd[itracks]));
          }
          //! for pions
          if (PDG_code_mc[id_from_mc_mpd[itracks]]==111)
          {
            h_pt_mpd_pi->Fill(TMath::Abs(signed_pt_mpd[itracks]));
          }
      		
   		   }

      }
      //! the same for mc-simulations  
      for (Long64_t itracks=0;itracks<n_tracks_mc;itracks++)
      {
        if(pt_mc[itracks]<=3 && TMath::Abs(eta_mc[itracks]) < 1.5  && mother_ID_mc[itracks] == -1) 
       	{
          if(PDG_code_mc[itracks]==2212)
          {
            h_pt_mc_p->Fill(pt_mc[itracks]);
          } 

          if(PDG_code_mc[itracks]==130)
          {
            h_pt_mc_k->Fill(pt_mc[itracks]);
          }

          if(PDG_code_mc[itracks]==111)
          {
            h_pt_mc_p->Fill(pt_mc[itracks]);
          } 

   		  }

      }

    }
   h_eff_p->Divide(h_pt_mpd_p, h_pt_mc_p);
   h_eff_k->Divide(h_pt_mpd_k, h_pt_mc_k);
   h_eff_pi->Divide(h_pt_mpd_pi, h_pt_mc_pi);

   //! for graphs
   h_eff_p->SetMarkerStyle(kDot);
   h_eff_k->SetMarkerStyle(kCircle);
   h_eff_pi->SetMarkerStyle(kStar);
   h_eff_p->SetMarkerColor(kRed);
   h_eff_k->SetMarkerColor(kGreen);
   h_eff_pi->SetMarkerColor(kBlue);

   TCanvas *c1 = new TCanvas("efficiency of p_T","efficiency of p_T", 600, 500);
   TH1F *f = (TH1F*)c1->DrawFrame(0,-1,3,7);
   f->GetXaxis()->SetTitle("p_T [GeV]");
   f->GetYaxis()->SetTitle("efficiency");
   h_eff_p->Draw();
   h_eff_k->Draw("same");
   h_eff_pi->Draw("same");

   TLegend *legend = new TLegend(0.55,0.65,0.16,0.82);
   legend->SetHeader("Efficiencies of particle reconstruction");
   legend->AddEntry(h_eff_p,"protons","lep");
   legend->AddEntry(h_eff_k,"kaons","lep");
   legend->AddEntry(h_eff_pi,"pions","lep");
   legend->Draw();

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
   output->Close();

}
