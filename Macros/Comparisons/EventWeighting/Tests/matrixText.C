void matrixText(){

  Double_t xVals[] = {1.0, 2.0, 3.0};
  Double_t yVals[] = {5.0, 6.0, 7.0};
  /*Double_t zVals[] = {1.0, 2.0, 3.0,
		      4.0, 5.0, 6.0,
		      7.0, 8.0, 9.0};//*/


  Double_t zVals[] = {3.0, 0.0, 0.0,
		      0.0, 2.0, 0.0,
		      0.0, 0.0, 1.0};
  
  TVectorD x; x.Use(3, xVals);
  TVectorD y; y.Use(3, yVals);
  TVectorD l = x-y;
  
  TMatrixDSym z(3); z.SetMatrixArray(zVals);
  //z.Invert();
  TVectorD k = z*y;
  //TVectorD k; k.Mult(z,y);

  //cout << ROOT::Math::Dot(x, y) << endl;
  cout << x*k << " " << x*y << " " << l*x << endl;//*/
}
