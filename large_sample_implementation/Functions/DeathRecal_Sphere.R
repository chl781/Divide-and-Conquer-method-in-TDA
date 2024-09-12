# This aims at recovering the death time for sphere case

DeathRecal_Sphere<- function(p1,p2,p3,p4){
    a = p1[1] - p2[1]
    b = p1[2] - p2[2]
    c = p1[3] - p2[3]
    
    a1 = p3[1] - p4[1]
    b1 = p3[2] - p4[2]
    c1 = p3[3] - p3[3]
    a2 = p2[1] - p3[1]
    b2 = p2[2] - p3[2]
    c2 = p2[3] - p3[3]
    
    A = p1[1] * p1[1] - p2[1] * p2[1]
    B = p1[2] * p1[2] - p2[2] * p2[2]
    C = p1[3] * p1[3] - p2[3] * p2[3]
    A1 = p3[1] * p3[1] - p4[1] * p4[1]
    B1 = p3[2] * p3[2] - p4[2] * p4[2]
    C1 = p3[3] * p3[3] - p4[3] * p4[3]
    A2 = p2[1] * p2[1] - p3[1] * p3[1]
    B2 = p2[2] * p2[2] - p3[2] * p3[2]
    C2 = p2[3] * p2[3] - p3[3] * p3[3]
    P = (A + B + C) /2
    Q = (A1 + B1 + C1) / 2
    R = (A2 + B2 + C2) / 2
    
    # D是系数行列式，利用克拉默法则
    D = a*b1*c2 + a2*b*c1 + c*a1*b2 - (a2*b1*c + a1*b*c2 + a*b2*c1)
    Dx = P*b1*c2 + b*c1*R + c*Q*b2 - (c*b1*R + P*c1*b2 + Q*b*c2)
    Dy = a*Q*c2 + P*c1*a2 + c*a1*R - (c*Q*a2 + a*c1*R + c2*P*a1)
    Dz = a*b1*R + b*Q*a2 + P*a1*b2 - (a2*b1*P + a*Q*b2 + R*b*a1)
    
    if(D == 0){
      print("4 points on the same plane" )
      return (-1)
    }else{
      centre=c(0,0,0)
      centre[1]=(Dx/D)
      centre[2]=(Dy/D)
      centre[3]=(Dz/D)
      radius = sqrt((p1[1]-centre[1])*(p1[1]-centre[1]) +
                      (p1[2]-centre[2])*(p1[2]-centre[2]) +
                      (p1[3]-centre[3])*(p1[3]-centre[3]))
      return(radius*sqrt(8/3))
    }
}
