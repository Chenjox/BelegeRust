use crate::schrauben::{ISOSchraube, self};




///
/// 
/// 
/// ```
///       width
///      |----|
///     
///      +----+  -
///      |    |  |
///    - |  P |  | height
///    | |    |  |
/// P.y| |    |  |
///    - o----+  -
///     
///      |--|
///      P.x
/// ```
pub struct Stahlblech {
  dicke: f64,  // in mm^2
  height: f64, // in mm^2
  width: f64,  // in mm^2
  festigkeit: f64 // in `N/mm^2`
}

impl Stahlblech {
  pub fn new(dicke: f64, height: f64, width: f64, festigkeit: f64) -> Self {
    return Stahlblech {
      dicke,
      height,
      width,
      festigkeit
    };
  }

  pub fn schwerpunkt(&self) -> [f64;2] {
    return [self.width/2., self.height/2.];
  }
}

pub struct Schraubenverbindung {
  blech: Stahlblech,
  schrauben: Vec<([f64;2],ISOSchraube)>
}

pub type Schraubengruppe = Vec<([f64;2],ISOSchraube)>;

impl Schraubenverbindung {
    pub fn new(blech: Stahlblech, schrauben: Schraubengruppe) -> Self {
      return Schraubenverbindung {
        blech,
        schrauben
      };
    }
    /// Mindestschraubenabstände p1, e1, p2, e2
    pub fn randabstaende(&self,idx: usize) -> [f64; 4] {
      let (coords, schraub) = &self.schrauben[idx];
      let mut p1 = f64::MAX;
      let mut p2 = f64::MAX;
      for cor in &self.schrauben {
        let cor = cor.0;
        if &cor == coords {
          continue;
        }
        p1 = p1.min((coords[0] - cor[0]).abs());
        p2 = p2.min((coords[1] - cor[1]).abs());
      }
      let e1 = ((self.blech.width - coords[0] ).abs()).min(coords[0]);
      let e2 = ((self.blech.height - coords[1]).abs()).min(coords[1]);

      return [p1,e1,p2,e2];
    }
}


