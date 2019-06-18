

#include<vector>
#include<list>

#include "plotcollection.h"
#include "complexMatrix.h"
#include "TFile.h"
#include "fitResult.h"
#include "TTree.h"
#include "TTreeCache.h"
#include <iostream>




int main(){


	std::cout << "start" << std::endl;
	double a;
	std::cin >> a;

	rpwa::fitResult resultIn;
	std::cout << sizeof(resultIn) << std::endl;


//	{
//		std::vector<std::string> files = {
//				"fits_tPrime-0.600-1.000.root", "fits_tPrime-0.140-0.210.root",
//				"fits_tPrime-0.210-0.600.root", "fits_tPrime-1.000-1.300.root",
//				"fits_tPrime-0.100-0.140.root" };
//
//		files = { "fits_tPrime-0.100-0.160.root", "fits_tPrime-0.160-0.260.root" };
//
////		std::vector<rpwa::fitResult> results;
//		std::list<rpwa::fitResult> results;
////		results.reserve(100000);
//
//		rpwa::fitResult resultIn;
//		rpwa::fitResult* resultInPtr = &resultIn;
//		for (const auto& fname : files) {
//			TFile* f = new TFile(fname.c_str());
//			TTree* t = (TTree*) f->Get("pwa");
//			t->SetBranchAddress("fitResult_v2", &resultInPtr);
//
////			results.reserve(results.size() + t->GetEntries());
//			for (int i = 0; i < t->GetEntries(); ++i) {
//				t->GetEntry(i);
//				results.push_back(resultIn);
//			}
////			TTreeCache* c = dynamic_cast<TTreeCache*>(f->GetCacheRead(t));
////			c->ResetCache();
////			c->Print();
//			t->PrintCacheStats();
//			f->Close();
//		}
////		std::cout << results.size() <<  '\t' << results.capacity()  << std::endl;
//		std::cout << "wait" << std::endl;
//		std::cin >> a;
//	}


//	{
//		std::vector<rpwa::complexMatrix> ms;
//		for(int i = 0; i < 10000; ++i){
//			rpwa::complexMatrix m(100,100);
//			for(int i = 0; i < 100; ++i){
//				for(int j = 0; j < 100; ++j){
//					m(i,j) = i*j;
//				}
//			}
//			ms.push_back(m);
//		}
//	std::cout << "wait" << std::endl;
//	std::cin >> a;
//	}


//	{
//		TFile* f = new TFile("fits_tPrime-0.100-0.700.root");
//		TTree* t = (TTree*) f->Get("pwa");
//		rpwa::fitResult resultIn;
//		rpwa::fitResult* resultInPtr = &resultIn;
//		t->SetBranchAddress("fitResult_v2", &resultInPtr);
//		t->GetEntry(0);
//		f->Close();
//		resultIn.stripCovMatrix();
//		resultIn.stripIntegralMatrices();
//		std::vector<rpwa::fitResult> rs;
//		rs.reserve(100000);
//		for(int i = 0; i < 100000; ++i){
//			rpwa::fitResult r;
//			r.fill(resultIn);
//			rs.push_back(r);
//		}
//	std::cout << "wait" << std::endl;
//	std::cin >> a;
//	}

//	{
//		std::vector<TMatrixT<Double_t> > test;
//		TMatrixT<Double_t> m(100, 100);
//		for(int i = 0; i < 1000; ++i){
//			TMatrixT<Double_t> n(100, 100);
//			n = m;
//			m[0][0] = i;
//			test.push_back(n);
//		}
//	std::cout << "wait" << std::endl;
//	std::cin >> a;
//	}

	{
//		typedef rpwa::complexMatrix Mat;
		typedef boost::numeric::ublas::matrix<std::complex<double>, boost::numeric::ublas::row_major  > Mat; //!
		std::vector<Mat > test;
		for(int i = 0; i < 1000; ++i){
			Mat m(0, 0);
			test.push_back(m);
		}
	std::cout << "wait" << std::endl;
	std::cin >> a;
	for(auto& m: test){ m.resize(100,100); for(int i=0; i<100; ++i){ for( int j=0; j < 100; ++j) m(i,j) = i*j;}; }
	std::cout << "wait" << std::endl;
	std::cin >> a;
	test.reserve(test.size()*2);
	std::cout << "wait" << std::endl;
	std::cin >> a;
	}

//	{
//		std::vector<double* > test;
//		for(int i = 0; i < 1000; ++i){
//			double* d = new double[100*100];
////			m.resizeTo(0,0);
//			test.push_back(d);
//		}
//	std::cout << "wait" << std::endl;
//	std::cin >> a;
//	test.reserve(test.size()*2);
//	std::cout << "wait" << std::endl;
//	std::cin >> a;
//	}



//	{
//		TFile* f = new TFile("out.root", "RECREATE");
//		TMatrixT<Double_t> m(100, 100);
//		TTree* t = new TTree("test", "test");
//		t->Branch("mat", &m);
//		for(int i = 0; i < 100000; ++i){
//			t->Fill();
//			m[i%100][i%100] = i*i;
//		}
//		f->Write();
//		f->Close();
//	std::cout << "wait" << std::endl;
//	std::cin >> a;
//	}

//	{
//		std::vector<std::string> files = {"out.root", "outII.root"};
//		for (const auto& fn: files){
//			TFile* f = new TFile(fn.c_str(), "READ");
//			TMatrixT<Double_t> m(100, 100);
//			TMatrixT<Double_t>* mPtr = &m;
//			TTree* t = (TTree*)f->Get("test");
//			t->SetBranchAddress("mat", &mPtr);
//			for(int i = 0; i < t->GetEntries(); ++i){
//				t->GetEntry(i);
//			}
//			f->Close();
//			std::cout << "wait" << std::endl;
//			std::cin >> a;
//		}
//	}

	std::cout << "done" << std::endl;
	std::cin >> a;

}
