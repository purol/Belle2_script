void CL_calculator(){

	const double CL = 0.95;
	const double step = 0.01;
	double Area = 0;
	double xmin = 7.66004;
	double xmax = 7.66004;

	TF1* f3 = new TF1("f3", "exp(-0.5*((x-7.66004)/19.4381)^2)", 0, 1000);

	const double total_Area = f3->Integral(0, 1000);

	while (total_Area*CL > Area) {
		xmin = xmin - step;
		xmax = xmax + step;
		if (xmin < 0) xmin = 0;
		Area = f3->Integral(xmin, xmax);
	}
	printf("xmin: %lf\nxmax: %lf\n", xmin, xmax);

	TCanvas* c1 = new TCanvas("c1", "Graph Draw Options", 200, 10, 600, 400);

	TF1* f1 = new TF1("f1", "exp(-0.5*((x-7.66004)/19.4381)^2)", 0, 100);
	f1->SetLineColor(1);
	f1->SetTitle("");
	f1->Draw();

	TF1* f2 = new TF1("f2", "exp(-0.5*((x-7.66004)/19.4381)^2)", xmin, xmax);
	f2->SetLineColor(1);
	f2->SetFillColor(15);
	f2->SetFillStyle(3004);
	f2->Draw("SAME");
}
