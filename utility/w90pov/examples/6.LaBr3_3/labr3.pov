// std povray files
#version 3.7;
global_settings {
assumed_gamma 1.0
max_trace_level 25
ambient_light rgb<1, 1, 1>*1
photons {
spacing 0.005
autostop 0
jitter 0
}
} 

#include "colors.inc"
#include "math.inc"
#include "glass.inc"
plane { <0, 0, 1>, -0.0
     pigment {color rgb 2.0 } 
    //pigment {
     // checker color Red, color Blue
    //}

photons {
  refraction on
  reflection on
  }
}

 sky_sphere {
    pigment {
      gradient y
      color_map {
        [ 0.5  color CornflowerBlue ]
        [ 1.0  color MidnightBlue ]
      }
      scale 2
      translate -1
    }
  }
 

// wannier stuff
#include "mydefs.inc"
#include "unitcell.inc"
#include "densities.inc"
#include "blobs.inc"
