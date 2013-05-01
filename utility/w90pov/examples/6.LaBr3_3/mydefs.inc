// degree of interpolation
#declare interp= 2;

 // Bravais lattice
 #declare a1 = <     3.9824000,    -6.8977191,     0.0000000>;
 #declare a2 = <     3.9824000,     6.8977191,     0.0000000>;
 #declare a3 = <     0.0000000,     0.0000000,     4.5119000>;
#declare lookpos=<3.9824,   2.2992,   1.1280>;
#declare Width = 400;
#declare Height = 400;
#declare minScreenDimension = 400;
#declare Scale =    6;
#declare Ratio = Scale * Width / Height;
#declare campos=vnormalize(a1)*Scale;
#declare RIGHT=<Ratio,0,0>;
#declare UP=vnormalize(vcross(campos,RIGHT))*Scale;
camera{
 //orthographic
 location campos+(a1+a2+a3)/2.0
 right RIGHT
 up UP
 sky   UP
 angle 30
 look_at lookpos
 }
 background { color rgb<1.0,1.0,1.0>}

 light_source {  lookpos+a3+0.3*a1 rgb <1.0,1.0,1.0> 
  photons {
  reflection on
  refraction on
}}


//***********************************************
// macros for common shapes
//***********************************************

#macro CW_angle (COLOR,A)
   #local RGBFT = color COLOR;
   #local R = (RGBFT.red);
   #local G = (RGBFT.green);
   #local B = (RGBFT.blue);
   #local Min = min(R,min(G,B));
   #local Max = max(R,max(G,B));
   #local Span = Max-Min;
   #local H = CRGB2H (<R,G,B>, Max, Span);
   #local S = 0; #if (Max!=0) #local S = Span/Max; #end

   #local P = <H+A,S,Max,(RGBFT.filter),(RGBFT.transmit)> ;

   #local HSVFT = color P ;
#local H = (HSVFT.red);
   #local S = (HSVFT.green);
   #local V = (HSVFT.blue);
   #local SatRGB = CH2RGB(H);
   #local RGB = ( ((1-S)*<1,1,1> + S*SatRGB) * V );
   rgb <RGB.red,RGB.green,RGB.blue,(HSVFT.filter),
       (HSVFT.transmit)>
#end

#default { finish {
 ambient .2 diffuse .6 specular 1 roughness .001 metallic}}

#macro atom(X,Y,Z,RADIUS,R,G,B)
 sphere{<X,Y,Z>,RADIUS
  pigment{rgb<R,G,B> } finish { phong 0.7 phong_size 90 }}
#end

#macro bond1(X1,Y1,Z1,X2,Y2,Z2,RADIUS,R,G,B)
 cylinder{<X1,Y1,Z1>,<X2,Y2,Z2>,RADIUS
  pigment{rgb<R,G,B>}}
  sphere{<X1,Y1,Z1>,RADIUS
   pigment{rgb<R,G,B>}}
  sphere{<X2,Y2,Z2>,RADIUS
   pigment{rgb<R,G,B>}}
#end
#macro bond2p(p1,p2,RADIUS,col1,col2)
 #declare rc=(p1+p2)/2.0;
 cylinder{p1, rc, RADIUS
  pigment{rgb col1} finish { phong 0.7 phong_size 90 } }
 cylinder{rc, p2, RADIUS
  pigment{rgb col2}finish { phong 0.7 phong_size 90 }}
  sphere{p1,RADIUS
   pigment{rgb col1}finish { phong 0.7 phong_size 90 }}
  sphere{p2,RADIUS
   pigment{rgb col2}finish { phong 0.7 phong_size 90 }}
#end



// make cage
#declare P = 2.4116;
#declare thick = 0.05;
#declare R = 0.5;
#declare G = 0.5;
#declare B = 0.5;


// My blob
#macro elblob1(level, thenorm, fun, c1, c2, c3, corg, col, trans)
object{

 isosurface { 
    function { level/abs(thenorm) - fun( x, y, z)  }
    accuracy 0.00001
    contained_by { box { <0,0,0>,<1.0,1.0,1.0> } }
    max_gradient 90.000
    pigment { rgbt <col.red,col.green,col.blue,trans>}
    finish { phong 0.7 phong_size 90 }
    }
    matrix <      
     c1.x,c1.y,c1.z,
     c2.x,c2.y,c2.z,
     c3.x,c3.y,c3.z,
     0,0,0>
     #translate corg
}
#end


#macro elblob2(level, thenorm, fun, c1, c2, c3, corg, col, trans)
object{

 isosurface { 
    function { level/abs(thenorm) - fun( x, y, z)  }
    accuracy 0.00001
    contained_by { box { <0,0,0>,<1.0,1.0,1.0> } }
    max_gradient 90.000
    pigment { rgbf <col.red,col.green,col.blue,0.8>}

    //finish { phong 0.7 phong_size 90 }
    //pigment{ rgbf<0.98,0.98,0.98,0.8>}
    finish { diffuse 0.01
             reflection 0.2
             specular 0.8
             roughness 0.0003
             phong 0.7
             phong_size 200}
interior {
I_Glass_Exp(2)

fade_color Col_Red_03 }
photons{
  target
  reflection on
  refraction on
  collect off 
}

    }
    matrix <      
     c1.x,c1.y,c1.z,
     c2.x,c2.y,c2.z,
     c3.x,c3.y,c3.z,
     0,0,0>
     #translate corg
}
#end


// A pm blob
#macro elblobpm(level, maxp, minp, funp, funm, c1, c2, c3, corg, col, trans)
 elblob2(level, maxp, funp, c1, c2, c3, corg, col, trans)
 elblob2(level, minp, funm, c1, c2, c3, corg, CW_angle (col,60), trans)
#end
 
