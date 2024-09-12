# DeathRecal aiming at circle case. This is derived from formula on the Internet.

DeathRecal_Circle<- function(X,i,j,k){
  A=X[i,]
  B=X[j,]
  C=X[k,]
  x1 = A[1];
  y1 = A[2];
  x2 = B[1];
  y2 = B[2];
  x3 = C[1];
  y3 = C[2];
  e = 2 * (x2 - x1);
  f = 2 * (y2 - y1);
  g = x2*x2 - x1*x1 + y2*y2 - y1*y1;
  a = 2 * (x3 - x2);
  b = 2 * (y3 - y2);
  c = x3*x3 - x2*x2 + y3*y3 - y2*y2;
  
  X = (g*b - c*f) / (e*b - a*f);
  Y = (a*g - c*e) / (a*f - b*e);
  R = sqrt((X-x1)*(X-x1)+(Y-y1)*(Y-y1));
  return(sqrt(3)*R)
}